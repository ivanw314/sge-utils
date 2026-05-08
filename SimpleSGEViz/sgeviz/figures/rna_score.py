import re

import altair as alt
import pandas as pd


_PALETTE = [
    "#006616",  # Synonymous
    "#81B4C7",  # Missense
    "#ffcd3a",  # Stop Gained
    "#93C47D",  # UTR Variant
    "#888888",  # Stop Lost
    "#000000",  # Start Lost
    "#CFCFCF",  # Splice Region
]

_CONSEQUENCE_DOMAIN = [
    "Synonymous",
    "Missense",
    "Stop Gained",
    "UTR Variant",
    "Stop Lost",
    "Start Lost",
    "Splice Region",
]


def _natsort_key(s: str) -> list:
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r"(\d+)", s)]


def _prepare_df(scores_df: pd.DataFrame) -> pd.DataFrame:
    """Filter to SNVs with RNA scores and recode consequence labels."""
    df = scores_df.loc[
        (scores_df["var_type"] == "snv") & scores_df["RNA_score"].notna()
    ].copy()

    rename = {
        "score": "snv_score",
        "RNA_score": "RNA/DNA",
        "CDS_position": "CDSpos",
        "amino_acid_change": "AAsub",
        "consequence": "Consequence",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})

    df.loc[df["Consequence"].str.contains("missense", na=False), "Consequence"] = "Missense"
    df.loc[df["Consequence"] == "synonymous_variant", "Consequence"] = "Synonymous"
    df.loc[df["Consequence"] == "intron_variant", "Consequence"] = "Intron"
    df.loc[df["Consequence"] == "stop_gained", "Consequence"] = "Stop Gained"
    df.loc[df["Consequence"] == "stop_lost", "Consequence"] = "Stop Lost"
    df.loc[df["Consequence"].str.contains("site", na=False), "Consequence"] = "Canonical Splice"
    df.loc[df["Consequence"].str.contains("ing_var", na=False), "Consequence"] = "Splice Region"
    df.loc[df["Consequence"].str.contains("UTR", na=False), "Consequence"] = "UTR Variant"
    df.loc[df["Consequence"] == "start_lost", "Consequence"] = "Start Lost"

    return df


def make_scatter(
    scores_df: pd.DataFrame,
    thresholds: list,
    rna_threshold: float | None = None,
    gene: str = "",
) -> alt.Chart:
    """RNA score vs. fitness score scatter plot.

    Always generated when RNA_score is present. The RNA threshold vertical line
    is omitted when rna_threshold is None.
    """
    alt.data_transformers.disable_max_rows()
    df = _prepare_df(scores_df)

    tooltips = [alt.Tooltip("target:N", title="SGE Region")]
    if "AAsub" in df.columns:
        tooltips.append(alt.Tooltip("AAsub:N", title="Amino Acid Substitution"))

    scatter = (
        alt.Chart(df)
        .mark_point()
        .encode(
            x=alt.X(
                "RNA/DNA:Q",
                title="RNA Score",
                axis=alt.Axis(
                    labelFontSize=16,
                    titleFontSize=20,
                    values=[-8, -6, -4, -2, 0, 2],
                ),
                scale=alt.Scale(domain=[-8, 3]),
            ),
            y=alt.Y(
                "snv_score:Q",
                title="Fitness Score",
                axis=alt.Axis(
                    labelFontSize=16,
                    titleFontSize=20,
                    values=[-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3],
                ),
                scale=alt.Scale(domain=[-0.6, 0.3]),
            ),
            color=alt.Color(
                "Consequence:N",
                scale=alt.Scale(range=_PALETTE, domain=_CONSEQUENCE_DOMAIN),
                legend=alt.Legend(
                    titleFontSize=16,
                    labelFontSize=14,
                    orient="bottom",
                    columns=3,
                    direction="horizontal",
                    titleOrient="left",
                    titleAnchor="start",
                ),
            ),
            tooltip=tooltips,
        )
        .properties(
            width=400,
            height=300,
            title=alt.TitleParams(
                text=f"Fitness Score vs. RNA Score{' (' + gene + ')' if gene else ''}",
                fontSize=22,
            ),
        )
    )

    nf_line = (
        alt.Chart(pd.DataFrame({"v": [thresholds[0]]}))
        .mark_rule(color="black", strokeDash=[8, 8], strokeWidth=2)
        .encode(y="v:Q")
    )
    func_line = (
        alt.Chart(pd.DataFrame({"v": [thresholds[1]]}))
        .mark_rule(color="#888888", strokeDash=[8, 8], strokeWidth=2)
        .encode(y="v:Q")
    )

    layers = [scatter, nf_line, func_line]

    if rna_threshold is not None:
        rna_line = (
            alt.Chart(pd.DataFrame({"v": [rna_threshold]}))
            .mark_rule(color="red", strokeWidth=2)
            .encode(x="v:Q")
        )
        layers.append(rna_line)

    return (
        alt.layer(*layers)
        .configure_axis(grid=False)
        .configure_view(stroke=None)
        .interactive()
    )


def make_stem_plot(
    scores_df: pd.DataFrame,
    rna_threshold: float,
    gene: str = "",
) -> alt.Chart:
    """Faceted stem plot of RNA scores per exon.

    Stems connect the RNA threshold line to any variant that is both
    functionally abnormal and RNA-low. Requires RNA_consequence column.
    """
    alt.data_transformers.disable_max_rows()
    df = _prepare_df(scores_df)

    if "CDSpos" in df.columns:
        df["CDSpos"] = pd.to_numeric(df["CDSpos"], errors="coerce")
        df = df.dropna(subset=["CDSpos"])

    df = df.loc[df["RNA/DNA"] >= -9].copy()

    consequence_sort = ["Functionally Normal", "Functionally Abnormal", "Indeterminate"]
    df.loc[df["functional_consequence"] == "functionally_abnormal", "functional_consequence"] = "Functionally Abnormal"
    df.loc[df["functional_consequence"] == "functionally_normal", "functional_consequence"] = "Functionally Normal"
    df.loc[df["functional_consequence"] == "indeterminate", "functional_consequence"] = "Indeterminate"

    df["min_threshold"] = rna_threshold
    df["rule_value"] = float("nan")
    stem_mask = (df["functional_consequence"] == "Functionally Abnormal") & (df["RNA_consequence"] == "low")
    df.loc[stem_mask, "rule_value"] = df.loc[stem_mask, "RNA/DNA"]

    exons = sorted(df["exon"].dropna().unique().tolist(), key=_natsort_key)

    tooltips = [
        alt.Tooltip("target:N", title="SGE Target"),
        alt.Tooltip("Consequence:N", title="Consequence"),
        alt.Tooltip("snv_score:Q", title="Fitness Score"),
    ]
    if "AAsub" in df.columns:
        tooltips.insert(1, alt.Tooltip("AAsub:N", title="Amino Acid Sub"))

    base = alt.Chart(df)

    dots = base.mark_point(filled=True, size=75).encode(
        x=alt.X(
            "CDSpos:Q",
            title="CDS Position",
            axis=alt.Axis(labelFontSize=20, titleFontSize=22),
            scale=alt.Scale(zero=False, padding=1),
        ),
        y=alt.Y(
            "RNA/DNA:Q",
            title="RNA Score",
            axis=alt.Axis(
                labelFontSize=20,
                titleFontSize=22,
                values=[-8, -6, -4, -2, 0, 2],
            ),
            scale=alt.Scale(domain=[-8, 3], padding=5),
        ),
        color=alt.Color(
            "Consequence:N",
            scale=alt.Scale(range=_PALETTE, domain=_CONSEQUENCE_DOMAIN),
            legend=alt.Legend(
                titleFontSize=22,
                labelFontSize=20,
                orient="top",
                columns=3,
                direction="horizontal",
                titleOrient="top",
                titleAnchor="start",
            ),
        ),
        shape=alt.Shape(
            "functional_consequence:N",
            sort=consequence_sort,
            legend=alt.Legend(
                title="Functional Consequence",
                titleLimit=500,
                titleFontSize=22,
                labelFontSize=20,
                labelLimit=500,
                orient="top",
                columns=1,
                direction="horizontal",
                titleOrient="top",
                titleAnchor="start",
            ),
        ),
        tooltip=tooltips,
    ).properties(width=600, height=150)

    stems = base.mark_rule().encode(
        x=alt.X("CDSpos:Q", scale=alt.Scale(zero=False, padding=1)),
        y=alt.Y("min_threshold:Q"),
        y2=alt.Y2("rule_value:Q"),
    )

    rna_line = (
        alt.Chart(pd.DataFrame({"y": [rna_threshold]}))
        .mark_rule(color="red")
        .encode(y=alt.Y("y:Q"))
    )

    return (
        (stems + dots + rna_line)
        .facet(facet=alt.Facet("exon:N", sort=exons), columns=2)
        .resolve_scale(x="independent", y="independent")
        .resolve_axis(x="independent", y="independent")
        .configure_axis(grid=False)
        .configure_view(stroke=None)
        .configure_header(labelFontSize=20, labelFontWeight="bold", title=None)
    )
