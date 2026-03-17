import altair as alt
import pandas as pd


def make_plot(df: pd.DataFrame, gene: str = "") -> alt.Chart:
    """Generate a 2D binned heatmap of log10(Allele Frequency) vs. SGE fitness score.

    Each bin is colored by the count of variants it contains. If both gnomAD and
    Regeneron datasets are present, the plot is faceted by dataset.

    Args:
        df: Combined dataframe from process.load_allele_freqs, containing
            'score', 'log_AF', and 'dataset' columns.
        gene: Gene name used in the plot title.
    """
    alt.data_transformers.disable_max_rows()

    n = len(df)
    title = f"Allele Frequency vs. Fitness Score{' (' + gene + ')' if gene else ''} (n = {n})"

    chart = alt.Chart(df).mark_rect().encode(
        x=alt.X(
            "score:Q",
            bin=alt.Bin(maxbins=50),
            axis=alt.Axis(
                title="Fitness Score",
                titleFontSize=20,
                labelFontSize=18,
                values=[-0.4, -0.2, 0],
            ),
        ),
        y=alt.Y(
            "log_AF:Q",
            bin=alt.Bin(maxbins=20),
            axis=alt.Axis(
                title="log10(Allele Frequency)",
                titleFontSize=20,
                labelFontSize=18,
                values=[0, -1, -2, -3, -4, -5, -6],
            ),
        ),
        color=alt.Color(
            "count():Q",
            scale=alt.Scale(scheme="lighttealblue", domain=[1, 400]),
            legend=alt.Legend(
                title="# of Variants", titleFontSize=16, labelFontSize=14
            ),
        ),
    ).properties(
        width=300,
        height=250,
        title=alt.TitleParams(text=title, fontSize=22),
    ).configure_axis(
        grid=False
    ).configure_view(
        stroke=None
    )

    return chart
