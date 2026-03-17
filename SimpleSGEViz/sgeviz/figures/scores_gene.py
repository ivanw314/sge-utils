import altair as alt
import pandas as pd
from natsort import natsorted

from .base import PALETTE, VARIANT_TYPES


def make_plot(df: pd.DataFrame, thresholds: list, gene: str = "") -> alt.Chart:
    """Generate a faceted scatter plot of SGE scores across gene exons.

    Each exon gets its own panel with an independent x-axis (genomic coordinates).
    Threshold lines mark the non-functional and functional score boundaries.
    """
    alt.data_transformers.disable_max_rows()

    n = len(df)
    df = df.copy()
    df["exon"] = df["exon"].transform(lambda x: "Exon " + x.split("_")[1][1:])
    df.loc[df["var_type"] == "snv", "var_type"] = "SNV"
    df.loc[df["var_type"] == "3bp_del", "var_type"] = "Deletion"

    exon_order = natsorted(df["exon"].unique().tolist())
    chromosome = df["chrom"].iloc[0] if "chrom" in df.columns else ""

    chart = alt.Chart(df).mark_point().encode(
        x=alt.X(
            "pos:Q",
            title=f"Genomic Coordinate ({chromosome})",
            axis=alt.Axis(titleFontSize=20, labelFontSize=18),
            scale=alt.Scale(zero=False),
        ),
        y=alt.Y(
            "score:Q",
            title="Fitness Score",
            axis=alt.Axis(titleFontSize=20, labelFontSize=18),
        ),
        color=alt.Color(
            "Consequence:N",
            scale=alt.Scale(domain=VARIANT_TYPES, range=PALETTE),
            legend=alt.Legend(
                symbolStrokeWidth=2,
                titleFontSize=20,
                labelFontSize=18,
                orient="bottom",
                columns=3,
                direction="horizontal",
            ),
        ),
        shape=alt.Shape(
            "var_type:N",
            legend=alt.Legend(
                title="Variant Type",
                symbolStrokeColor="black",
                symbolStrokeWidth=2,
                symbolFillColor="black",
                titleFontSize=20,
                labelFontSize=18,
                orient="bottom",
            ),
        ),
    ).properties(width=600, height=150)

    nf_line = alt.Chart(df).mark_rule(
        color="black", strokeDash=[8, 8], strokeWidth=1
    ).encode(y=alt.datum(thresholds[0]))

    func_line = alt.Chart(df).mark_rule(
        color="#888888", strokeDash=[8, 8], strokeWidth=1
    ).encode(y=alt.datum(thresholds[1]))

    return (
        alt.layer(chart, nf_line, func_line)
        .facet(
            facet=alt.Facet(
                "exon:N",
                sort=exon_order,
                title=f"{gene + ' ' if gene else ''}(n = {n})",
            ),
            columns=2,
        )
        .resolve_scale(x="independent")
        .configure_header(
            labelFontSize=18,
            titleFontSize=20,
            labelFontWeight="bold",
        )
    )
