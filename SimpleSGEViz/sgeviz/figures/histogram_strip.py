import math

import pandas as pd
import altair as alt

from .base import PALETTE, VARIANT_TYPES


def _score_tick_values(scores: pd.Series) -> list:
    """X-axis tick positions spanning the actual score range.

    Tries steps [0.5, 0.2, 0.1, 0.05, 0.02] from coarsest to finest and
    picks the coarsest one that still produces at least 6 ticks, so narrow
    datasets get sensible density rather than 2-3 widely-spaced marks.
    """
    data_min, data_max = scores.min(), scores.max()
    for step in [0.5, 0.2, 0.1, 0.05, 0.02]:
        lo = math.floor(round(data_min / step, 6)) * step
        hi = math.ceil(round(data_max / step, 6)) * step
        n = round((hi - lo) / step) + 1
        if n >= 6:
            break
    return [round(lo + i * step, 10) for i in range(n)]


def _count_tick_values(scores: pd.Series, nbins: int = 50) -> list:
    """Y-axis tick positions for the histogram based on approximate peak bin count."""
    max_count = int(pd.cut(scores, bins=nbins).value_counts().max())
    for step in [50, 100, 200, 500, 1000, 2000, 5000]:
        if math.ceil(max_count / step) <= 8:
            break
    hi = math.ceil(max_count / step) * step
    return list(range(0, hi + step, step))


def _threshold_rules(thresholds: list):
    """Return (non-functional rule, functional rule) Altair mark_rule charts."""
    nf_line = alt.Chart(
        pd.DataFrame({"Fitness Score": [thresholds[0]]})
    ).mark_rule(color="black", strokeDash=[8, 8], strokeWidth=2).encode(
        x="Fitness Score:Q"
    )

    func_line = alt.Chart(
        pd.DataFrame({"Fitness Score": [thresholds[1]]})
    ).mark_rule(color="#888888", strokeDash=[8, 8], strokeWidth=2).encode(
        x="Fitness Score:Q"
    )

    return nf_line, func_line


def make_figures(df: pd.DataFrame, thresholds: list, gene: str = ""):
    """Build histogram and strip plot charts.

    Returns:
        histogram: Altair layer chart
        stripplot: Altair layer chart
    """
    alt.data_transformers.disable_max_rows()

    n = len(df)
    x_ticks = _score_tick_values(df["score"])
    y_ticks = _count_tick_values(df["score"])
    nf_line, func_line = _threshold_rules(thresholds)
    selection = alt.selection_point(fields=["Consequence"], bind="legend")
    zoom_hist = alt.selection_interval(bind="scales", name="zoom_hist")
    zoom_strip = alt.selection_interval(bind="scales", name="zoom_strip")

    histogram = alt.Chart(df).mark_bar().encode(
        alt.X(
            "score",
            bin=alt.Bin(maxbins=50),
            axis=alt.Axis(
                title="",
                labelFontSize=16,
                titleFontSize=20,
                values=x_ticks,
            ),
        ),
        alt.Y(
            "count()",
            axis=alt.Axis(
                title="Number of Variants",
                labelFontSize=16,
                titleFontSize=20,
                values=y_ticks,
            ),
        ),
        color=alt.Color(
            "Consequence:N",
            scale=alt.Scale(range=PALETTE, domain=VARIANT_TYPES),
            legend=alt.Legend(
                titleFontSize=16, labelFontSize=14,
                orient="none", legendX=10, legendY=5,
                columns=2,
            ),
        ),
        opacity=alt.condition(selection, alt.value(1), alt.value(0.2)),
    ).add_params(selection, zoom_hist).properties(
        width=800,
        height=200,
        title=alt.TitleParams(
            text=f"Distribution of SGE Scores{' (' + gene + ')' if gene else ''} (n = {n})",
            fontSize=22,
        ),
    )

    histogram = alt.layer(histogram, nf_line, func_line).resolve_scale(y="shared")

    df = df.copy()
    df["Consequence"] = pd.Categorical(
        df["Consequence"], categories=VARIANT_TYPES, ordered=True
    )

    stripplot = alt.Chart(df).mark_tick(opacity=1).encode(
        x=alt.X(
            "score:Q",
            axis=alt.Axis(
                title="Fitness Score",
                values=x_ticks,
                titleFontSize=20,
                labelFontSize=16,
            ),
        ),
        y=alt.Y(
            "Consequence:N",
            sort=VARIANT_TYPES,
            axis=alt.Axis(title="", labelFontSize=16),
        ),
        color=alt.Color(
            "Consequence:N",
            legend=None,
            scale=alt.Scale(range=PALETTE, domain=VARIANT_TYPES),
        ),
    ).add_params(zoom_strip).properties(width=800, height=200)

    stripplot = alt.layer(stripplot, nf_line, func_line).resolve_scale(y="shared")

    return histogram, stripplot


def combine(histogram: alt.Chart, stripplot: alt.Chart) -> alt.Chart:
    """Stack histogram and strip plot vertically with shared x-axis and clean styling."""
    return (
        (histogram & stripplot)
        .configure_axis(grid=False)
        .configure_view(stroke=None)
        .resolve_scale(x="shared")
    )
