from pathlib import Path

import altair as alt
import pandas as pd

_DOMAIN_PALETTE = [
    "#B9DBF4",  # light blue
    "#C8DBC8",  # light green
    "#FF9A00",  # orange
    "#ffc976",  # amber
    "#D5A6BD",  # pink
    "#A4C2F4",  # periwinkle
    "#FFD966",  # gold
    "#B6D7A8",  # sage green
]

_DEL_TYPES = ["Inframe Indel", "Stop Gained"]
_DEL_PALETTE = ["black", "#ffc007"]
_DEL_CONSEQUENCE_MAP = {
    # process.py already maps inframe_indel → "Inframe Indel" and stop_gained → "Stop Gained";
    # merge the rarer cases into "Inframe Indel" for the del panel display.
    "Start Lost": "Inframe Indel",
    "Stop Lost": "Inframe Indel",
}

_AA_ORDER = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
    "*", "Mis. Min.", "Mis. Mean",
]

_VEP_COLS = {
    "am_score": "AlphaMissense",
    "revel_score": "REVEL",
    "MutPred2": "MutPred2",
}


def _load_domains(path, prot_length: int) -> pd.DataFrame:
    """Load a domain annotation file and return a single flat segment DataFrame.

    Expected columns: region_name, aa_residues (format "start-end").
    Optional columns:
      color (hex)  – auto-assigned from _DOMAIN_PALETTE if absent.
      tier (int)   – 0 = main domain (default), 1 = sub-feature.

    Tier-1 sub-features are inserted into the tier-0 strip in place: any tier-0
    segment overlapping a sub-feature is split into [before, sub-feature, after],
    so both colors appear inline in a single strip.

    Returns a DataFrame with start, end, label, color columns.
    """
    suffix = Path(path).suffix.lower()
    raw = pd.read_excel(path) if suffix in (".xlsx", ".xls") else pd.read_csv(path)

    raw[["start", "end"]] = raw["aa_residues"].str.split("-", expand=True).astype(int)
    raw = raw.sort_values("start").reset_index(drop=True)

    if "tier" not in raw.columns:
        raw["tier"] = 0

    if "color" not in raw.columns:
        palette = _DOMAIN_PALETTE * ((len(raw) // len(_DOMAIN_PALETTE)) + 1)
        raw["color"] = palette[: len(raw)]

    # Build tier-0 base segments with gray gap fills
    tier0 = raw.loc[raw["tier"] == 0].sort_values("start").reset_index(drop=True)
    segments = []
    pos = 1
    for _, row in tier0.iterrows():
        if pos < row["start"]:
            segments.append({"start": pos, "end": row["start"], "label": "", "color": "#CCCCCC", "tier": 0})
        segments.append({
            "start": row["start"],
            "end": row["end"] + 1,
            "label": row["region_name"],
            "color": row["color"],
            "tier": 0,
        })
        pos = row["end"] + 1
    if pos <= prot_length:
        segments.append({"start": pos, "end": prot_length + 1, "label": "", "color": "#CCCCCC", "tier": 0})

    # Insert tier-1 sub-features by splitting any overlapping parent segment
    tier1 = raw.loc[raw["tier"] == 1].sort_values("start").reset_index(drop=True)
    for _, sub in tier1.iterrows():
        sub_end = sub["end"] + 1  # convert inclusive end to exclusive
        new_segs = []
        for seg in segments:
            if sub["start"] >= seg["end"] or sub_end <= seg["start"]:
                new_segs.append(seg)
                continue
            # Split: piece before sub-feature
            if seg["start"] < sub["start"]:
                new_segs.append({"start": seg["start"], "end": sub["start"],
                                 "label": seg["label"], "color": seg["color"], "tier": seg["tier"]})
            # The sub-feature itself (clipped to parent bounds)
            new_segs.append({
                "start": max(seg["start"], sub["start"]),
                "end": min(seg["end"], sub_end),
                "label": sub["region_name"],
                "color": sub["color"],
                "tier": 1,
            })
            # Piece after sub-feature
            if seg["end"] > sub_end:
                new_segs.append({"start": sub_end, "end": seg["end"],
                                 "label": seg["label"], "color": seg["color"], "tier": seg["tier"]})
        segments = new_segs

    return pd.DataFrame(segments)



def _prep_aa_exon_df(aa_exon_df: "pd.DataFrame | None") -> "pd.DataFrame | None":
    """Convert aa_start/aa_end columns to the start/end/label/center format
    expected by _make_exon_cartoon. Returns None if input is None."""
    if aa_exon_df is None or aa_exon_df.empty:
        return None
    df = aa_exon_df.rename(columns={"aa_start": "start", "aa_end": "end"}).copy()
    if "exon_num" in df.columns:
        df["label"] = df["exon_num"].astype(str)
    else:
        df["label"] = [str(i + 1) for i in range(len(df))]
    df["center"] = (df["start"] + df["end"]) / 2
    return df[["start", "end", "label", "center"]]


def _make_exon_cartoon(exon_df: pd.DataFrame, prot_length: int, width: int) -> alt.Chart:
    """White rect strip showing exon boundaries with sequential exon numbers."""
    RECT_H = 26
    x_scale = alt.Scale(domain=[1, prot_length + 1])

    rects = (
        alt.Chart(exon_df)
        .mark_rect(stroke="black", strokeWidth=1.5, color="white")
        .encode(
            x=alt.X("start:Q", axis=None, scale=x_scale),
            x2="end:Q",
            tooltip=[
                alt.Tooltip("label:N", title="Exon"),
                alt.Tooltip("start:Q", title="AA Start"),
                alt.Tooltip("end:Q", title="AA End"),
            ],
        )
        .properties(width=width, height=RECT_H)
    )

    labels = (
        alt.Chart(exon_df)
        .mark_text(fontSize=14, fontWeight="bold", baseline="middle", align="center")
        .encode(
            x=alt.X("center:Q", axis=None, scale=x_scale),
            text="label:N",
        )
        .properties(width=width, height=RECT_H)
    )

    return alt.layer(rects, labels)


def _make_domain_cartoon(segments_df: pd.DataFrame, prot_length: int, width: int) -> alt.Chart:
    """Colored rect strip with labels in a zone above the strip.

    Labels that are too wide for their segment, or that would overlap a
    neighbour, are suppressed. The segment is still rendered and has a tooltip.
    """
    LABEL_H = 22   # pixels reserved above the colored strip for labels
    RECT_H = 28    # height of the colored strip itself
    TOTAL_H = LABEL_H + RECT_H
    CHAR_W = 9     # approximate px per character at 16pt bold
    PAD = 4        # horizontal padding on each side of a label
    x_scale = alt.Scale(domain=[1, prot_length + 1])

    df = segments_df.copy()
    df["center"] = (df["start"] + df["end"]) / 2

    label_df = df.loc[df["label"] != ""].copy()
    # y position: center of the label zone above the strip
    label_df["y_mid"] = 1.0 - (LABEL_H / 2) / TOTAL_H

    # Pixel geometry for overlap detection
    label_df["center_px"] = (label_df["center"] - 1) / prot_length * width
    label_df["seg_px"] = (label_df["end"] - label_df["start"]) / prot_length * width
    label_df["txt_px"] = label_df["label"].str.len() * CHAR_W + PAD * 2

    # Drop labels that don't fit within their own segment (tier-1 labels are exempt)
    is_tier1 = label_df["tier"] == 1 if "tier" in label_df.columns else pd.Series(False, index=label_df.index)
    label_df = label_df.loc[(label_df["txt_px"] <= label_df["seg_px"]) | is_tier1].reset_index(drop=True)

    # Tier-priority label suppression:
    #   1. Keep all tier-1 labels that fit (they take priority).
    #   2. Greedy left-to-right pass for tier-0 labels, skipping any that
    #      would overlap a reserved tier-1 label slot.
    if not label_df.empty:
        has_tier = "tier" in label_df.columns

        # Collect reserved pixel ranges from tier-1 labels
        tier1_ranges = []
        if has_tier:
            for _, row in label_df.loc[label_df["tier"] == 1].iterrows():
                tier1_ranges.append((
                    row["center_px"] - row["txt_px"] / 2,
                    row["center_px"] + row["txt_px"] / 2,
                ))

        def _overlaps_tier1(left_px, right_px):
            return any(left_px < r1 and right_px > l1 for l1, r1 in tier1_ranges)

        keep = []
        # Always keep tier-1 labels
        if has_tier:
            keep.extend(label_df.loc[label_df["tier"] == 1].index.tolist())

        # Greedy pass for tier-0 labels
        tier0_labels = label_df.loc[label_df["tier"] == 0] if has_tier else label_df
        last_right_px = -float("inf")
        for _, row in tier0_labels.sort_values("center_px").iterrows():
            left_px = row["center_px"] - row["txt_px"] / 2
            right_px = row["center_px"] + row["txt_px"] / 2
            if left_px >= last_right_px and not _overlaps_tier1(left_px, right_px):
                keep.append(row.name)
                last_right_px = right_px

        label_df = label_df.loc[label_df.index.isin(keep)]

    rects = (
        alt.Chart(df)
        .mark_rect(stroke="black", strokeWidth=0.5)
        .encode(
            x=alt.X("start:Q", axis=None, scale=x_scale),
            x2="end:Q",
            y=alt.value(LABEL_H),
            y2=alt.value(TOTAL_H),
            color=alt.Color("color:N", scale=None, legend=None),
            tooltip=[
                alt.Tooltip("label:N", title="Domain"),
                alt.Tooltip("start:Q", title="Start"),
                alt.Tooltip("end:Q", title="End"),
            ],
        )
        .properties(width=width, height=TOTAL_H)
    )

    labels = (
        alt.Chart(label_df)
        .mark_text(fontSize=16, fontWeight="bold", baseline="middle", align="center")
        .encode(
            x=alt.X("center:Q", axis=None, scale=x_scale),
            y=alt.Y("y_mid:Q", scale=alt.Scale(domain=[0, 1]), axis=None),
            text="label:N",
        )
        .properties(width=width, height=TOTAL_H)
    )

    return alt.layer(rects, labels).resolve_scale(y="independent")


def _make_del_panel(
    del_df: pd.DataFrame, thresholds, prot_length: int, width: int
) -> alt.Chart:
    """Scatter plot of 3bp deletion fitness scores by protein position."""
    del_df = del_df.copy()
    if "amino_acid_change" in del_df.columns:
        del_df = del_df.loc[~del_df["amino_acid_change"].isin(["---"])]

    del_df["_cds_start"] = del_df["CDS_position"].str.split("-").str[0]
    del_df["_cds_end"] = del_df["CDS_position"].str.split("-").str[1]
    del_df = del_df.loc[
        ~del_df["_cds_start"].isin(["?"]) & ~del_df["_cds_end"].isin(["?"])
    ]
    del_df["_cds_start"] = del_df["_cds_start"].astype(int)
    del_df["ps_aa_start"] = ((del_df["_cds_start"] + 2) / 3).round(2)
    del_df = del_df.loc[del_df["ps_aa_start"] <= prot_length]
    del_df["Consequence"] = del_df["Consequence"].replace(_DEL_CONSEQUENCE_MAP)

    y_min = min(-0.5, del_df["score"].min())
    y_max = max(0.1, del_df["score"].max())

    scatter = (
        alt.Chart(del_df)
        .mark_point(strokeWidth=3, size=75, opacity=1)
        .encode(
            x=alt.X(
                "ps_aa_start:Q",
                title="",
                axis=alt.Axis(ticks=False, labels=False),
                scale=alt.Scale(domain=[1, prot_length + 1], nice=False),
            ),
            y=alt.Y(
                "score:Q",
                title="Fitness Score",
                axis=alt.Axis(titleFontSize=20, labelFontSize=18),
                scale=alt.Scale(domain=[y_min, y_max]),
            ),
            color=alt.Color(
                "Consequence:N",
                scale=alt.Scale(domain=_DEL_TYPES, range=_DEL_PALETTE),
                legend=alt.Legend(
                    title="Consequence", labelFontSize=18, titleFontSize=20, orient="right"
                ),
            ),
            shape=alt.Shape(
                "Consequence:N",
                scale=alt.Scale(domain=_DEL_TYPES, range=["square", "triangle"]),
                legend=None,
            ),
            order=alt.Order("Consequence:N", sort="ascending"),
        )
        .properties(width=width, height=150)
    )

    layers = [scatter]
    if thresholds is not None:
        layers.append(
            alt.Chart(del_df)
            .mark_rule(color="red", strokeDash=[8, 8], strokeWidth=1)
            .encode(y=alt.datum(thresholds[0]))
        )
        layers.append(
            alt.Chart(del_df)
            .mark_rule(color="blue", strokeDash=[8, 8], strokeWidth=1)
            .encode(y=alt.datum(thresholds[1]))
        )
    return alt.layer(*layers)


def make_plot(
    df: pd.DataFrame,
    gene: str = "",
    thresholds=None,
    domains_path=None,
    protein_length: int | None = None,
    px_per_aa: int = 3,
    aa_exon_df: "pd.DataFrame | None" = None,
) -> alt.Chart:
    """Generate amino acid substitution heatmap of SGE fitness scores.

    X-axis: amino acid position (per-residue bins).
    Y-axis: amino acid substitution, sorted by 1-letter code then *, Mis. Min., Mis. Mean.
    Color: fitness score (bluepurple reversed, clamped to [-0.2, 0]).

    If am_score, revel_score, or MutPred2 columns are present in df, a VEP predictor
    sub-panel is appended below the main heatmap with an independent color scale.

    Args:
        df: scores dataframe from process.load_scores (must contain amino_acid_change).
        gene: Gene name used in the plot title.
        protein_length: Known full length of the protein. If provided and larger than
            the maximum position observed in the data, it will be used instead (the
            dataset may be incomplete and not cover all residues).
        px_per_aa: Pixels allocated per amino acid column (default 3). Reduce to
            produce a narrower figure.
    """
    alt.data_transformers.disable_max_rows()

    snv_df = df.loc[df["var_type"] == "snv"].copy()
    snv_df = snv_df.loc[~snv_df["amino_acid_change"].isin(["-", "---"])]
    snv_df = snv_df.dropna(subset=["amino_acid_change"])

    # Parse amino acid change string: e.g. 'A123G' -> og_AA='A', AA_change='G', AApos=123
    snv_df["og_AA"] = snv_df["amino_acid_change"].str[0]
    snv_df["AA_change"] = snv_df["amino_acid_change"].str[-1]
    snv_df["AApos"] = pd.to_numeric(
        snv_df["amino_acid_change"].str[1:-1], errors="coerce"
    )
    snv_df = snv_df.dropna(subset=["AApos"])
    snv_df["AApos"] = snv_df["AApos"].astype(int)
    snv_df = snv_df.loc[snv_df["og_AA"] != "*"]

    if "max_SpliceAI" in snv_df.columns:
        snv_df = snv_df.loc[snv_df["max_SpliceAI"] <= 0.2]

    prot_length = snv_df["AApos"].max()
    if protein_length is not None and protein_length > prot_length:
        prot_length = protein_length
    n = len(snv_df)
    n_del = int((df["var_type"] == "3bp_del").sum()) if "var_type" in df.columns else 0
    width = px_per_aa * prot_length
    height_per_row = 25

    # Missense min/mean rows (exclude stop gained)
    mis_df = snv_df.loc[snv_df["Consequence"] != "Stop Gained"]
    min_df = mis_df.groupby("AApos")["score"].min().reset_index()
    min_df["AA_change"] = "Mis. Min."
    mean_df = mis_df.groupby("AApos")["score"].mean().reset_index()
    mean_df["AA_change"] = "Mis. Mean"

    plot_df = pd.concat(
        [snv_df[["AApos", "AA_change", "score"]], min_df, mean_df],
        ignore_index=True,
    )

    gene_label = f" ({gene})" if gene else ""
    count_label = f"n = {n_del} deletions, {n} SNVs" if n_del > 0 else f"n = {n}"
    title = f"Deletion and Heatmap{gene_label} ({count_label})"

    vep_cols_present = {k: v for k, v in _VEP_COLS.items() if k in snv_df.columns}
    has_vep = bool(vep_cols_present)

    heatmap = alt.Chart(plot_df).mark_rect().encode(
        x=alt.X(
            "AApos:Q",
            title="" if has_vep else "Amino Acid Position",
            axis=alt.Axis(
                labels=not has_vep,
                ticks=not has_vep,
                values=[] if has_vep else list(range(0, prot_length, 25)),
                titleFontSize=22,
                labelFontSize=18,
            ),
            scale=alt.Scale(domain=[1, prot_length + 1]),
            bin=alt.Bin(maxbins=prot_length + 1, minstep=1),
        ),
        y=alt.Y(
            "AA_change:N",
            title="Amino Acid Substitution",
            sort=_AA_ORDER,
            axis=alt.Axis(labelFontSize=18, titleFontSize=20),
        ),
        color=alt.Color(
            "score:Q",
            title="Fitness Score",
            scale=alt.Scale(domain=[-0.2, 0], clamp=True, scheme="bluepurple", reverse=True),
            legend=alt.Legend(titleFontSize=18, labelFontSize=16),
        ),
    ).properties(width=width, height=height_per_row * len(_AA_ORDER))

    if has_vep:
        vep_summary = (
            snv_df.groupby("AApos")[list(vep_cols_present.keys())]
            .mean()
            .reset_index()
            .rename(columns=vep_cols_present)
        )
        vep_df = vep_summary.melt(id_vars=["AApos"], var_name="Predictor", value_name="score")
        vep_order = list(vep_cols_present.values())

        vep_map = alt.Chart(vep_df).mark_rect().encode(
            x=alt.X(
                "AApos:Q",
                title="Amino Acid Position",
                axis=alt.Axis(
                    values=list(range(0, prot_length, 25)),
                    titleFontSize=22,
                    labelFontSize=18,
                ),
                scale=alt.Scale(domain=[1, prot_length + 1]),
                bin=alt.Bin(maxbins=prot_length + 1, minstep=1),
            ),
            y=alt.Y(
                "Predictor:N",
                title="",
                sort=vep_order,
                axis=alt.Axis(labelFontSize=16),
            ),
            color=alt.Color(
                "score:Q",
                title="Predictor Score",
                scale=alt.Scale(domain=[0, 1], clamp=True, scheme="bluepurple"),
                legend=alt.Legend(
                    titleFontSize=18,
                    labelFontSize=16,
                    gradientLength=height_per_row * len(vep_cols_present),
                    values=[0, 1.0],
                ),
            ),
        ).properties(width=width, height=height_per_row * len(vep_cols_present))

        lower = alt.vconcat(heatmap, vep_map, spacing=9).resolve_scale(color="independent")
    else:
        lower = heatmap

    # Clip exon aa_end to prot_length + 1 so the stop codon doesn't expand the shared x scale
    if aa_exon_df is not None:
        aa_exon_df = aa_exon_df.copy()
        aa_exon_df["aa_end"] = aa_exon_df["aa_end"].clip(upper=prot_length + 1)

    # Build panel stack top-to-bottom: domains → deletions → heatmap(+vep)
    panels = []

    if domains_path is not None:
        segments_df = _load_domains(domains_path, prot_length)
        domain_chart = _make_domain_cartoon(segments_df, prot_length, width)
        exon_df = _prep_aa_exon_df(aa_exon_df)
        if exon_df is not None:
            exon_chart = _make_exon_cartoon(exon_df, prot_length, width)
            panels.append(alt.vconcat(domain_chart, exon_chart, spacing=0))
        else:
            panels.append(domain_chart)
    elif aa_exon_df is not None:
        panels.append(_make_exon_cartoon(_prep_aa_exon_df(aa_exon_df), prot_length, width))

    has_del = (
        "var_type" in df.columns
        and "CDS_position" in df.columns
        and (df["var_type"] == "3bp_del").any()
    )
    if has_del:
        del_df = df.loc[df["var_type"] == "3bp_del"].copy()
        panels.append(_make_del_panel(del_df, thresholds, prot_length, width))

    panels.append(lower)

    if len(panels) == 1:
        chart = panels[0]
    else:
        chart = alt.vconcat(*panels, spacing=9).resolve_scale(
            x="shared", color="independent", shape="independent"
        )

    return (
        chart
        .properties(title=alt.TitleParams(text=title, fontSize=22, anchor="middle"))
        .configure_axis(grid=False)
        .configure_view(stroke=None)
    )
