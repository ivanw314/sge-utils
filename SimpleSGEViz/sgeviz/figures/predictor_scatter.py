"""Predictor score vs. SGE fitness score scatter plots.

Generates one scatter panel per available predictor (AlphaMissense, REVEL,
CADD, MutPred2), arranged two-per-row.  Each panel shows SNVs colored by
variant consequence, with vertical lines at the SGE thresholds and horizontal
lines at published benign/pathogenic thresholds for missense predictors.
"""

import altair as alt
import pandas as pd

from .base import PALETTE, VARIANT_TYPES

# Column name → display label, in the order panels will appear
_PREDICTORS = [
    ("am_score",    "AlphaMissense"),
    ("revel_score", "REVEL"),
    ("cadd_score",  "CADD"),
    ("MutPred2",    "MutPred2"),
]

# Missense-specific predictors: filtered to max_SpliceAI <= 0.2 when available
_MISSENSE_ONLY = {"am_score", "revel_score", "MutPred2"}

# Published moderate-evidence thresholds: (benign_max, pathogenic_min)
# Source: Bergquist et al. 2025, Nat Genet (doi:10.1038/s41588-024-01908-2)
_PRED_THRESHOLDS = {
    "am_score":    (0.099, 0.906),
    "revel_score": (0.183, 0.773),
    "MutPred2":    (0.197, 0.829),
}

_PANEL_W = 350
_PANEL_H = 250


def make_plot(df: pd.DataFrame, thresholds: list, gene: str = "") -> alt.Chart | None:
    """Generate predictor vs. SGE fitness score scatter panels.

    One scatter panel is produced per available predictor, arranged two-per-row.
    Each panel shows:
      - SNV variants colored by standardized consequence (filtered by
        max_SpliceAI <= 0.2 for missense-specific predictors)
      - Vertical dashed lines at SGE non-functional (black) and functional
        (gray) thresholds
      - Horizontal dashed lines at published benign (blue) and pathogenic
        (red) thresholds for AlphaMissense, REVEL, and MutPred2

    Returns None if none of the expected predictor columns are present.
    """
    alt.data_transformers.disable_max_rows()

    snv_df = df.loc[df["var_type"] == "snv"].copy()

    panels = []
    for col, display in _PREDICTORS:
        if col not in snv_df.columns:
            continue

        panel_df = snv_df.copy()
        if col in _MISSENSE_ONLY and "max_SpliceAI" in panel_df.columns:
            panel_df = panel_df.loc[panel_df["max_SpliceAI"] <= 0.2]
        panel_df = panel_df.dropna(subset=[col])
        if panel_df.empty:
            continue

        # Rename predictor column to display name for clean axis label
        panel_df = panel_df.rename(columns={col: display})

        scatter = alt.Chart(panel_df).mark_circle(
            opacity=0.7, size=20
        ).encode(
            x=alt.X(
                "score:Q",
                axis=alt.Axis(
                    title="Fitness Score",
                    titleFontSize=16,
                    labelFontSize=14,
                ),
                scale=alt.Scale(zero=False),
            ),
            y=alt.Y(
                f"{display}:Q",
                axis=alt.Axis(
                    title=display,
                    titleFontSize=16,
                    labelFontSize=14,
                ),
            ),
            color=alt.Color(
                "Consequence:N",
                scale=alt.Scale(domain=VARIANT_TYPES, range=PALETTE),
                legend=alt.Legend(
                    title="Consequence",
                    titleFontSize=14,
                    labelFontSize=12,
                    symbolOpacity=1,
                ),
            ),
            tooltip=[
                alt.Tooltip("pos_id:N", title="Position ID"),
                alt.Tooltip("Consequence:N"),
                alt.Tooltip("score:Q", title="Fitness Score", format=".3f"),
                alt.Tooltip(f"{display}:Q", title=display, format=".3f"),
            ],
        ).properties(
            width=_PANEL_W,
            height=_PANEL_H,
            title=alt.TitleParams(
                text=f"{gene + ' — ' if gene else ''}{display} vs. Fitness Score",
                fontSize=15,
            ),
        ).interactive(name=f"pan_{col}")

        # SGE threshold lines (vertical)
        nf_line = alt.Chart(panel_df).mark_rule(
            color="black", strokeDash=[8, 8], strokeWidth=1.5
        ).encode(x=alt.datum(thresholds[0]))

        func_line = alt.Chart(panel_df).mark_rule(
            color="#888888", strokeDash=[8, 8], strokeWidth=1.5
        ).encode(x=alt.datum(thresholds[1]))

        layers = [scatter, nf_line, func_line]

        # Predictor threshold lines (horizontal) where published thresholds exist
        if col in _PRED_THRESHOLDS:
            benign_thresh, path_thresh = _PRED_THRESHOLDS[col]
            layers.append(
                alt.Chart(panel_df).mark_rule(
                    color="#1a7abf", strokeDash=[4, 4], strokeWidth=1.5
                ).encode(y=alt.datum(benign_thresh))
            )
            layers.append(
                alt.Chart(panel_df).mark_rule(
                    color="#d62728", strokeDash=[4, 4], strokeWidth=1.5
                ).encode(y=alt.datum(path_thresh))
            )

        panels.append(alt.layer(*layers))

    if not panels:
        return None

    rows = [
        alt.hconcat(*panels[i:i + 2])
        for i in range(0, len(panels), 2)
    ]
    return (
        alt.vconcat(*rows)
        .configure_axis(grid=False)
        .configure_view(stroke=None)
    )
