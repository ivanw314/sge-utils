from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd


_DOMAIN_COLORS = [
    "#B9DBF4", "#C8DBC8", "#F6BF93", "#D5D0F2",
    "#018571", "#D35400", "#2980B9", "#C0392B",
]


def _build_vcoords(
    exon_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    intron_vw: float = 0.04,
) -> tuple:
    """Map genomic exon coordinates to a normalised visual [0, ~1] x-axis.

    CDS exon widths are proportional to their genomic bp. UTR sub-segments
    and introns each receive *intron_vw* of visual space.

    Returns (segments, gv) where segments is a list of dicts with keys
    kind/gstart/gend/vstart/vend/exon, and gv(genomic_pos) maps a genomic
    position to its visual x-coordinate.
    """
    meta   = dict(zip(meta_df["type"], meta_df["info"]))
    strand = str(meta.get("strand", "plus")).lower()
    atg    = int(meta["atg"])
    stop   = int(meta["stop"])

    orig = exon_df.sort_values("start").reset_index(drop=True)
    if strand == "minus":
        exons = orig.copy()
        exons["start"] = -orig["end"].values
        exons["end"]   = -orig["start"].values
        atg, stop = -atg, -stop
    else:
        exons = orig.copy()
    exons = exons.sort_values("start").reset_index(drop=True)

    def _split(gs, ge):
        subs, p = [], gs
        if p < atg:
            subs.append((p, min(ge, atg), "utr")); p = min(ge, atg)
        if p < ge and p < stop:
            subs.append((p, min(ge, stop), "cds")); p = min(ge, stop)
        if p < ge:
            subs.append((p, ge, "utr"))
        return subs

    all_subs = [
        (i, g0, g1, k, str(row.get("exon", i + 1)))
        for i, row in exons.iterrows()
        for g0, g1, k in _split(int(row["start"]), int(row["end"]))
    ]

    utr_max_vw    = 3 * intron_vw
    utr_bps       = [g1 - g0 for _, g0, g1, k, __ in all_subs if k == "utr"]
    utr_px_per_bp = utr_max_vw / max(utr_bps) if utr_bps else 0.0
    utr_vw        = lambda bp: max(utr_px_per_bp * bp, intron_vw)

    cds_bp     = sum(g1 - g0 for _, g0, g1, k, __ in all_subs if k == "cds")
    compressed = (len(exons) - 1) * intron_vw + sum(
        utr_vw(g1 - g0) for _, g0, g1, k, __ in all_subs if k == "utr"
    )
    cds_vw    = max(0.1, 1.0 - compressed)
    px_per_bp = cds_vw / max(cds_bp, 1)

    segments: list = []
    vpos = 0.0
    prev_idx = None

    for exon_idx, g0, g1, kind, label in all_subs:
        if prev_idx is not None and exon_idx != prev_idx:
            i0 = int(exons.iloc[prev_idx]["end"])
            i1 = int(exons.iloc[exon_idx]["start"])
            segments.append({
                "kind": "intron", "gstart": i0, "gend": i1,
                "vstart": vpos, "vend": vpos + intron_vw, "exon": None,
            })
            vpos += intron_vw

        vw = utr_vw(g1 - g0) if kind == "utr" else px_per_bp * (g1 - g0)
        segments.append({
            "kind": kind, "gstart": g0, "gend": g1,
            "vstart": vpos, "vend": vpos + vw, "exon": label,
        })
        vpos += vw
        prev_idx = exon_idx

    def gv(gpos: float) -> float:
        p = -gpos if strand == "minus" else gpos
        for seg in segments:
            if seg["gstart"] <= p <= seg["gend"]:
                span = seg["gend"] - seg["gstart"]
                t = (p - seg["gstart"]) / span if span else 0.0
                return seg["vstart"] + t * (seg["vend"] - seg["vstart"])
        return segments[0]["vstart"] if p < segments[0]["gstart"] else segments[-1]["vend"]

    return segments, gv


def _build_aa_to_genomic_map(exon_df: pd.DataFrame, meta_df: pd.DataFrame):
    """Return aa_to_gen(aa_pos) -> genomic coordinate for domain coloring."""
    meta   = dict(zip(meta_df["type"], meta_df["info"]))
    strand = str(meta.get("strand", "plus")).lower()
    atg    = int(meta["atg"])
    stop   = int(meta["stop"])
    exons  = exon_df.sort_values("start").reset_index(drop=True)

    if strand == "plus":
        cds_segs = [
            (max(int(r["start"]), atg), min(int(r["end"]), stop))
            for _, r in exons.iterrows()
            if min(int(r["end"]), stop) > max(int(r["start"]), atg)
        ]
        def aa_to_gen(aa_pos):
            target, acc = (aa_pos - 1) * 3, 0
            for gs, ge in cds_segs:
                if acc + (ge - gs) > target:
                    return gs + (target - acc)
                acc += ge - gs
            return cds_segs[-1][1]
    else:
        cds_segs = [
            (max(int(r["start"]), stop), min(int(r["end"]), atg))
            for _, r in exons.sort_values("start", ascending=False).iterrows()
            if min(int(r["end"]), atg) > max(int(r["start"]), stop)
        ]
        def aa_to_gen(aa_pos):
            target, acc = (aa_pos - 1) * 3, 0
            for cds_lo, cds_hi in cds_segs:
                if acc + (cds_hi - cds_lo) > target:
                    return cds_hi - (target - acc)
                acc += cds_hi - cds_lo
            return cds_segs[-1][0]

    return aa_to_gen


def _load_domains(domains_path: Path) -> pd.DataFrame:
    """Load pipeline domain file and return a DataFrame with start, end, domain columns."""
    df = pd.read_csv(domains_path) if domains_path.suffix == ".csv" else pd.read_excel(domains_path)
    if "tier" in df.columns:
        df = df[df["tier"] == 0]
    df[["start", "end"]] = df["aa_residues"].str.split("-", expand=True).astype(int)
    return df.rename(columns={"region_name": "domain"})[["start", "end", "domain"]]


def make_plot(
    gene: str,
    exon_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    lib_df: pd.DataFrame,
    scores_df: pd.DataFrame | None = None,
    domains_path: Path | None = None,
    exon_color: str = "#d0d0d0",
    hatch: str = "//////",
    fig_width: float = 10,
) -> plt.Figure:
    """Draw a single-track SGE library cartoon.

    CDS exon blocks are colored by protein domain where applicable (exon_color
    elsewhere). UTR sub-regions are drawn as thinner blocks. Introns are shown
    as backbone gaps. Library amplicon coverage is crosshatched over CDS/UTR
    regions. ATG and Stop markers are drawn above.

    When scores_df is provided, variant counts (SNVs, deletions, RNA
    measurements) are printed to the right of the gene name.
    """
    segments, gv = _build_vcoords(exon_df, meta_df)
    total_vw = segments[-1]["vend"]
    meta     = dict(zip(meta_df["type"], meta_df["info"]))

    domain_df = _load_domains(domains_path) if domains_path is not None else None
    has_dom   = domain_df is not None and not domain_df.empty

    BB          =  0.0
    CB,  CT     = -0.28,  0.28
    UB,  UT     = -0.17,  0.17
    LBL         = -0.40
    DOM_LBL     = -0.65
    MKR         =  0.40

    y_bot = DOM_LBL - 0.08 if has_dom else LBL - 0.08
    fig, ax = plt.subplots(figsize=(fig_width, 1.5 if has_dom else 1.2))

    # Backbone across intron gaps only
    for seg in segments:
        if seg["kind"] == "intron":
            ax.plot([seg["vstart"], seg["vend"]], [BB, BB], color="black", lw=1.2, zorder=1)

    # Base exon/UTR blocks and collect per-exon visual spans for labels
    exon_vspans: dict = {}
    for seg in segments:
        if seg["kind"] == "intron":
            continue
        name = seg["exon"]
        exon_vspans.setdefault(name, [seg["vstart"], seg["vend"]])
        exon_vspans[name][0] = min(exon_vspans[name][0], seg["vstart"])
        exon_vspans[name][1] = max(exon_vspans[name][1], seg["vend"])
        yb, yt = (CB, CT) if seg["kind"] == "cds" else (UB, UT)
        ax.add_patch(mpatches.Rectangle(
            (seg["vstart"], yb), seg["vend"] - seg["vstart"], yt - yb,
            facecolor=exon_color, edgecolor="black", lw=0.8, zorder=2,
        ))

    # Domain coloring overlay on CDS segments
    if has_dom:
        aa_to_gen = _build_aa_to_genomic_map(exon_df, meta_df)
        for i, (_, dom) in enumerate(domain_df.iterrows()):
            color = _DOMAIN_COLORS[i % len(_DOMAIN_COLORS)]
            v0 = gv(aa_to_gen(int(dom["start"])))
            v1 = gv(aa_to_gen(int(dom["end"])))
            vlo, vhi = min(v0, v1), max(v0, v1)
            dom_v_extents = []
            for seg in segments:
                if seg["kind"] != "cds":
                    continue
                ov0 = max(seg["vstart"], vlo)
                ov1 = min(seg["vend"],   vhi)
                if ov1 > ov0:
                    ax.add_patch(mpatches.Rectangle(
                        (ov0, CB), ov1 - ov0, CT - CB,
                        facecolor=color, edgecolor="black", lw=0.8, zorder=3,
                    ))
                    dom_v_extents.append((ov0, ov1))
            if dom_v_extents:
                label_x = (min(s[0] for s in dom_v_extents) + max(s[1] for s in dom_v_extents)) / 2
                ax.text(label_x, DOM_LBL, dom["domain"],
                        ha="center", va="top", fontsize=14, fontweight="bold", family="Arial")

    # Library amplicon crosshatch
    for _, row in lib_df.iterrows():
        v0 = gv(int(row["start"]))
        v1 = gv(int(row["end"]))
        vlo, vhi = min(v0, v1), max(v0, v1)
        for seg in segments:
            if seg["kind"] == "intron":
                continue
            ov0 = max(seg["vstart"], vlo)
            ov1 = min(seg["vend"],   vhi)
            if ov1 > ov0:
                yb, yt = (CB, CT) if seg["kind"] == "cds" else (UB, UT)
                ax.add_patch(mpatches.Rectangle(
                    (ov0, yb), ov1 - ov0, yt - yb,
                    facecolor="none", edgecolor="black", lw=0, hatch=hatch, zorder=4,
                ))

    # Exon number labels
    for name, (vs, ve) in sorted(exon_vspans.items(), key=lambda kv: kv[1][0]):
        lbl = name[1:] if name[:1] == "X" and name[1:].isdigit() else str(name)
        ax.text((vs + ve) / 2, LBL, lbl,
                ha="center", va="top", fontsize=14, fontweight="bold", family="Arial")

    # ATG / Stop markers
    for key, label in [("atg", "ATG"), ("stop", "Stop")]:
        if key in meta:
            vx = gv(int(meta[key]))
            ax.annotate("", xy=(vx, CT), xytext=(vx, MKR),
                        arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
            ax.text(vx, MKR + 0.04, label,
                    ha="center", va="bottom", fontsize=12, fontweight="bold", family="Arial")

    # Gene name
    ax.text(1.02, 0.90, gene, transform=ax.transAxes, clip_on=False,
            ha="left", va="center", fontsize=16, fontweight="bold", family="Arial")

    # Variant stats
    if scores_df is not None and not scores_df.empty:
        n_nts = scores_df["start"].nunique()
        n_snv = (scores_df["var_type"] == "snv").sum()
        n_del = (scores_df["var_type"] == "3bp_del").sum()
        lines = [
            f"{n_nts:,} nts targeted",
            f"{n_snv:,} SNVs",
            f"{n_del:,} 3-bp deletions",
        ]
        if "RNA_score" in scores_df.columns:
            n_rna = scores_df["RNA_score"].notna().sum()
            lines.append(f"{n_rna:,} RNA measurements")
        ax.text(1.02, 0.58, "\n".join(lines), transform=ax.transAxes, clip_on=False,
                ha="left", va="top", fontsize=14, family="Arial", linespacing=1.6)

    ax.set_xlim(-0.01, total_vw + 0.01)
    ax.set_ylim(y_bot, MKR + 0.25)
    ax.axis("off")
    plt.tight_layout()
    return fig
