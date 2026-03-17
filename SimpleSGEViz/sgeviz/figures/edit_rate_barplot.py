import re

import altair as alt
import pandas as pd


_REP_MAP = {
    "R1R4": "Rep. 1", "R1R2R3": "Rep. 1", "R1": "Rep. 1", "R4": "Rep. 1",
    "R2R5": "Rep. 2", "R4R5R6": "Rep. 2", "R2": "Rep. 2", "R5": "Rep. 2",
    "R3R6": "Rep. 3", "R7R8R9": "Rep. 3", "R3": "Rep. 3", "R6": "Rep. 3",
}

_REP_ORDER = ["Rep. 1", "Rep. 2", "Rep. 3"]


def _natsort_key(s: str) -> list:
    """Natural sort key: splits on digit runs so '10A' sorts after '9A'."""
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r"(\d+)", s)]


def make_plot(df: pd.DataFrame, gene: str = "") -> alt.Chart:
    """Return a faceted bar chart of library edit rates by target and replicate.

    Args:
        df: DataFrame with columns target_rep, edit_rate (raw TSV, not yet processed).
        gene: Gene name used in the figure title.

    Returns:
        A configured Altair Chart ready to save.
    """
    alt.data_transformers.disable_max_rows()

    df = df.copy()

    # Extract target exon label and replicate string from target_rep
    # e.g. CTCF_X10A_R1R4_D05 -> target='10A', rep='R1R4'
    df["target"] = df["target_rep"].transform(lambda x: x.split("X")[1].split("_")[0])
    df["rep"] = df["target_rep"].transform(lambda x: x.split("X")[1].split("_")[1])

    # Map replicate combos to display labels; fall back to raw string
    df["rep"] = df["rep"].map(_REP_MAP).fillna(df["rep"])

    sort_order = sorted(df["target"].unique().tolist(), key=_natsort_key)

    plot = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X(
                "rep:N",
                axis=alt.Axis(title="", labels=False, ticks=False),
            ),
            y=alt.Y(
                "edit_rate:Q",
                title="Lib. Edit Rate",
                axis=alt.Axis(labelFontSize=18, titleFontSize=20),
                scale=alt.Scale(domain=[0, 0.5]),
            ),
            column=alt.Column(
                "target:N",
                sort=sort_order,
                header=alt.Header(
                    title=f"SGE Target{' (' + gene + ')' if gene else ''}",
                    titleFontSize=22,
                    labelFontSize=20,
                ),
            ),
            color=alt.Color(
                "rep:N",
                sort=_REP_ORDER,
                legend=alt.Legend(
                    title="",
                    orient="bottom",
                    direction="horizontal",
                    labelFontSize=18,
                ),
            ),
        )
        .properties(width=35, height=200)
        .configure_facet(spacing=5)
        .configure_axis(grid=False)
    )

    return plot
