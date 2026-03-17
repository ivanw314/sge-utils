import altair as alt
import pandas as pd
from natsort import natsorted


_DAYS = ["D05", "D13", "D17"]
_REPS = ["R1", "R2", "R3"]


def _build_combos(columns: list) -> tuple[list, list]:
    """Return (rep_cols, combos) based on which day columns are present."""
    rep_cols = [
        f"{day} {rep}"
        for day in _DAYS
        for rep in _REPS
        if f"{day} {rep}" in columns
    ]
    combos = [
        (col1, col2)
        for i, col1 in enumerate(rep_cols)
        for col2 in rep_cols[i + 1 :]
        if col1[:3] == col2[:3]  # same day only
    ]
    return rep_cols, combos


def compute_correlations(counts_df: pd.DataFrame) -> pd.DataFrame:
    """Compute pairwise Pearson r between all replicate combinations per SGE target.

    Automatically detects which day timepoints (D05, D13, D17) are present and
    computes within-day pairwise correlations for each.

    Returns a long-form dataframe with columns: Targets, Tests, r_correlation.
    """
    rep_cols, combos = _build_combos(counts_df.columns.tolist())
    records = []

    for target, group in counts_df.groupby("target"):
        group_reps = group[rep_cols]
        for col1, col2 in combos:
            r = round(group_reps[col1].corr(group_reps[col2]), 3)
            records.append(
                {"Targets": target, "Tests": f"{col1} vs {col2}", "r_correlation": r}
            )

    return pd.DataFrame(records)


def make_heatmap(df: pd.DataFrame, gene: str = "") -> alt.Chart:
    """Generate Pearson r heatmap across SGE targets and replicate comparisons."""
    targets = natsorted(df["Targets"].unique().tolist())

    title = f"Correlation of Replicates{' (' + gene + ')' if gene else ''}"

    base = alt.Chart(
        df,
        title=alt.TitleParams(text=title, fontSize=32),
    ).encode(
        x=alt.X(
            "Tests:N",
            axis=alt.Axis(
                title="",
                titleFontSize=28,
                labelFontSize=24,
                labelLimit=300,
                labelAngle=45,
            ),
        ),
        y=alt.Y(
            "Targets:N",
            sort=targets,
            axis=alt.Axis(
                title="SGE Target Region", titleFontSize=28, labelFontSize=24
            ),
        ),
    )

    row_h = 34  # px per row: fits fontSize=20 text with ~14px breathing room
    height = len(targets) * row_h

    rect = base.mark_rect().encode(
        color=alt.Color(
            "r_correlation:Q",
            scale=alt.Scale(domain=[0.2, 1]),
            legend=alt.Legend(
                title="Pearson's r", titleFontSize=24, labelFontSize=22
            ),
        ),
        tooltip=[alt.Tooltip("r_correlation", title="Pearson's r: ")],
    ).properties(width=600, height=height)

    text_color = (
        alt.when(alt.datum.r_correlation > 0.5)
        .then(alt.value("white"))
        .otherwise(alt.value("black"))
    )

    text = base.mark_text(baseline="middle", fontSize=20).encode(
        text=alt.Text("r_correlation:Q", format="0.2f"),
        color=text_color,
    ).transform_filter("isValid(datum.r_correlation)")

    return rect + text
