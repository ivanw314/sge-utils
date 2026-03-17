import altair as alt
import numpy as np
import pandas as pd

from .base import PALETTE, VARIANT_TYPES


_CLASS_SORT = ["Benign", "Likely benign", "Uncertain significance", "Likely pathogenic", "Pathogenic"]


def _compute_roc(df: pd.DataFrame):
    """Compute ROC curve for SGE score classifying B/LB (negative) vs P/LP (positive).

    Lower SGE score → more non-functional → predicted pathogenic.
    Sweeps score thresholds ascending so pathogenic predictions accumulate first.

    Returns:
        roc_df: DataFrame with FPR and TPR columns, or None if insufficient data.
        auc_val: Area under the ROC curve (float), or None.
    """
    roc_df = df.loc[df["Germline classification"].isin(
        ["Benign", "Likely benign", "Pathogenic", "Likely pathogenic"]
    )].copy()

    roc_df["label"] = roc_df["Germline classification"].isin(
        ["Pathogenic", "Likely pathogenic"]
    ).astype(int)

    n_pos = int(roc_df["label"].sum())
    n_neg = len(roc_df) - n_pos
    if n_pos == 0 or n_neg == 0:
        return None, None

    roc_df = roc_df.sort_values("score").reset_index(drop=True)
    labels = roc_df["label"].values

    fprs, tprs = [0.0], [0.0]
    tp, fp = 0, 0
    for label in labels:
        if label == 1:
            tp += 1
        else:
            fp += 1
        tprs.append(tp / n_pos)
        fprs.append(fp / n_neg)

    auc_val = float(np.trapezoid(tprs, fprs))
    return pd.DataFrame({"FPR": fprs, "TPR": tprs}), auc_val


def make_strip(df: pd.DataFrame, thresholds: list, gene: str = "") -> alt.Chart:
    """Generate ClinVar classification strip plot.

    Per-variant SGE scores arranged by germline classification, colored by
    consequence type (same palette as other figures), with per-category counts
    annotated and threshold lines marking functional/non-functional boundaries.

    Args:
        df: merged ClinVar + scores dataframe from process.load_clinvar.
        thresholds: [non_functional_threshold, functional_threshold].
        gene: Gene name used in the plot title.
    """
    alt.data_transformers.disable_max_rows()

    n = len(df)
    title = f"{gene + ' ' if gene else ''}ClinVar Variants vs. SGE Score (n = {n})"

    nf_line = alt.Chart(
        pd.DataFrame({"x": [thresholds[0]]})
    ).mark_rule(color="black", strokeDash=[8, 8], strokeWidth=2).encode(x="x:Q")

    func_line = alt.Chart(
        pd.DataFrame({"x": [thresholds[1]]})
    ).mark_rule(color="#888888", strokeDash=[8, 8], strokeWidth=2).encode(x="x:Q")

    # ClinVar is SNV-only — exclude Inframe Indel from the legend
    _snv_types, _snv_palette = zip(
        *[(t, c) for t, c in zip(VARIANT_TYPES, PALETTE) if t != "Inframe Indel"]
    )

    strip = alt.Chart(df).mark_tick(opacity=1, thickness=2).encode(
        x=alt.X(
            "score:Q",
            title="Fitness Score",
            axis=alt.Axis(labelFontSize=16, titleFontSize=20),
        ),
        y=alt.Y(
            "Germline classification:N",
            title="",
            sort=_CLASS_SORT,
            axis=alt.Axis(labelFontSize=16, titleFontSize=20, labelLimit=300),
        ),
        color=alt.Color(
            "Consequence:N",
            scale=alt.Scale(domain=list(_snv_types), range=list(_snv_palette)),
            legend=alt.Legend(titleFontSize=16, labelFontSize=14),
        ),
        tooltip=["pos_id:N", "score:Q", "Consequence:N", "Germline classification:N"],
    ).properties(width=700, height=250)

    # Per-classification count annotations pinned to the left edge
    summary = (
        df.groupby("Germline classification")
        .size()
        .reset_index(name="count")
    )
    summary["text"] = "(n = " + summary["count"].astype(str) + ")"

    counts = alt.Chart(summary).mark_text(
        align="right", dx=-10, dy=20, fontSize=14, color="black"
    ).encode(
        y=alt.Y("Germline classification:N", sort=_CLASS_SORT),
        text="text:N",
        x=alt.value(0),
    )

    return (
        alt.layer(strip, counts, nf_line, func_line)
        .resolve_scale(y="shared")
        .properties(title=alt.TitleParams(text=title, fontSize=22))
        .configure_axis(grid=False)
        .configure_view(stroke=None)
    )


def make_roc(df: pd.DataFrame, gene: str = "") -> alt.Chart | None:
    """Generate ROC curve for SGE score classifying B/LB vs P/LP ClinVar variants.

    Returns None if there are insufficient variants in both classes.

    Args:
        df: merged ClinVar + scores dataframe from process.load_clinvar.
        gene: Gene name used in the plot title.
    """
    alt.data_transformers.disable_max_rows()

    roc_points, auc_val = _compute_roc(df)
    if roc_points is None:
        return None

    title = f"{gene + ' ' if gene else ''}ROC Curve: SGE Score (B/LB vs P/LP)"

    diagonal = pd.DataFrame({"FPR": [0.0, 1.0], "TPR": [0.0, 1.0]})
    diag_line = alt.Chart(diagonal).mark_line(
        color="gray", strokeDash=[4, 4], strokeWidth=1
    ).encode(x="FPR:Q", y="TPR:Q")

    roc_line = alt.Chart(roc_points).mark_line(
        color="orange", strokeWidth=2.5
    ).encode(
        x=alt.X(
            "FPR:Q",
            title="1 - Specificity (FPR)",
            scale=alt.Scale(domain=[0, 1]),
            axis=alt.Axis(labelFontSize=14, titleFontSize=16),
        ),
        y=alt.Y(
            "TPR:Q",
            title="Sensitivity (TPR)",
            scale=alt.Scale(domain=[0, 1]),
            axis=alt.Axis(labelFontSize=14, titleFontSize=16),
        ),
    )

    auc_label = alt.Chart(
        pd.DataFrame({"x": [0.55], "y": [0.08], "label": [f"AUC = {auc_val:.3f}"]})
    ).mark_text(fontSize=20, fontWeight="bold", align="left", color="black").encode(
        x="x:Q", y="y:Q", text="label:N"
    )

    return (
        alt.layer(diag_line, roc_line, auc_label)
        .properties(
            width=350,
            height=350,
            title=alt.TitleParams(text=title, fontSize=22),
        )
        .configure_axis(grid=False)
        .configure_view(stroke=None)
    )
