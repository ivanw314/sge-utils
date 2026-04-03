"""Genomic gene-structure and library-design cartoons.

The coordinate system uses a non-linear x-axis where:
  - CDS exon widths are proportional to their genomic length (equal px-per-bp).
  - Intron widths are linearly scaled between ``intron_min_vw`` and
    ``intron_max_vw`` pixels based on their genomic length, clamped to
    [``intron_min_bp``, ``intron_max_bp``].  This keeps very large introns
    from dominating the figure while still showing size differences between
    smaller ones.
  - UTR regions (5′ of ATG, 3′ of stop) are compressed the same way using
    ``utr_min_vw`` / ``utr_max_vw`` pixel limits.

A ``//`` marker is drawn on the backbone at the midpoint of every intron to
signal the compressed / non-linear scale.

Public functions
----------------
make_exon_cartoon   – exon structure track only.
make_library_cartoon – exon track stacked above a library-amplicon track.

Colors, ATG/stop positions, and strand are read from a metadata DataFrame
with columns ``['type', 'info']``.  Recognised ``type`` values:

  atg        genomic position of the start codon
  stop       genomic position of the stop codon
  strand     gene strand: 'plus' (default) or 'minus'
  exon_color hex color for exon rectangles  (default #2E86C1)
  lib_color  hex color for library amplicons (default #888888)

Minus-strand genes (``strand: minus``) are handled by negating all
coordinates internally so that the 5′ end of the transcript always appears
on the left.  All input coordinates (exon_df, lib_df, atg, stop) should be
provided in standard genomic convention (start < end, positive-strand
numbering).
"""

from __future__ import annotations

import altair as alt
import pandas as pd

_LEFT_MARGIN = 80   # px reserved on the left for row labels ("Exon", "Library")


# ── Coordinate transform helpers ──────────────────────────────────────────────

def _build_visual_segments(
    exon_df: pd.DataFrame,
    total_width: int,
    intron_min_bp: int,
    intron_max_bp: int,
    intron_min_vw: float,
    intron_max_vw: float,
    utr_min_bp: int,
    utr_max_bp: int,
    utr_min_vw: float,
    utr_max_vw: float,
    atg_pos: int,
    stop_pos: int,
) -> tuple[list[dict], float]:
    """Map genomic coordinates to visual x-coordinates.

    Each exon is split at the ATG and stop positions into UTR and CDS
    sub-segments.  UTR sub-segments receive a compressed visual width
    (scaled between *utr_min_vw* and *utr_max_vw*); CDS sub-segments share
    the remaining space proportionally.

    Returns ``(segments, total_vw)`` where each segment is a dict::

        {gstart, gend, vstart, vend, kind ('exon'|'utr'|'intron'), exon (str|None)}

    ``total_vw`` equals ``total_width`` by construction.
    """
    exons = exon_df.sort_values("start").reset_index(drop=True)
    n_exons = len(exons)

    def _scale(bp: int, min_bp: int, max_bp: int, min_vw: float, max_vw: float) -> float:
        t = max(0.0, min(1.0, (bp - min_bp) / max(max_bp - min_bp, 1)))
        return min_vw + t * (max_vw - min_vw)

    # Intron visual widths
    intron_vws: list[float] = []
    for idx in range(n_exons - 1):
        gap_bp = int(exons.iloc[idx + 1]["start"]) - int(exons.iloc[idx]["end"])
        intron_vws.append(_scale(gap_bp, intron_min_bp, intron_max_bp, intron_min_vw, intron_max_vw))

    # Split each exon into [5'UTR, CDS, 3'UTR] sub-segments
    def _split_exon(gstart: int, gend: int) -> list[tuple[int, int, str]]:
        subs: list[tuple[int, int, str]] = []
        pos = gstart
        if pos < atg_pos:
            end = min(gend, atg_pos)
            subs.append((pos, end, "utr"))
            pos = end
        if pos < gend and pos < stop_pos:
            end = min(gend, stop_pos)
            subs.append((pos, end, "exon"))
            pos = end
        if pos < gend:
            subs.append((pos, gend, "utr"))
        return subs

    # Pre-compute UTR visual widths and CDS bp total
    all_subs: list[tuple[int, int, int, str, str]] = []  # (exon_idx, g0, g1, kind, label)
    for idx in range(n_exons):
        row = exons.iloc[idx]
        for g0, g1, kind in _split_exon(int(row["start"]), int(row["end"])):
            all_subs.append((idx, g0, g1, kind, str(row["exon"])))

    utr_vw_map: dict[tuple[int, int, int], float] = {}
    total_cds_bp = 0
    for idx, g0, g1, kind, _ in all_subs:
        bp = g1 - g0
        if kind == "utr":
            utr_vw_map[(idx, g0, g1)] = _scale(bp, utr_min_bp, utr_max_bp, utr_min_vw, utr_max_vw)
        else:
            total_cds_bp += bp

    total_cds_vw = max(1.0, total_width - sum(intron_vws) - sum(utr_vw_map.values()))
    px_per_bp = total_cds_vw / max(total_cds_bp, 1)

    # Build segment list
    segments: list[dict] = []
    vpos = 0.0
    intron_idx = 0
    for idx in range(n_exons):
        row = exons.iloc[idx]
        for g0, g1, kind in _split_exon(int(row["start"]), int(row["end"])):
            bp = g1 - g0
            vw = utr_vw_map[(idx, g0, g1)] if kind == "utr" else px_per_bp * bp
            segments.append({
                "gstart": g0,
                "gend": g1,
                "vstart": vpos,
                "vend": vpos + vw,
                "kind": kind,
                "exon": str(row["exon"]),
            })
            vpos += vw
        if idx < n_exons - 1:
            segments.append({
                "gstart": int(row["end"]),
                "gend": int(exons.iloc[idx + 1]["start"]),
                "vstart": vpos,
                "vend": vpos + intron_vws[intron_idx],
                "kind": "intron",
                "exon": None,
            })
            vpos += intron_vws[intron_idx]
            intron_idx += 1

    return segments, vpos


def _gv(pos: float, segments: list[dict]) -> float:
    """Interpolate a genomic position to its visual x-coordinate."""
    for seg in segments:
        if seg["gstart"] <= pos <= seg["gend"]:
            span = seg["gend"] - seg["gstart"]
            if span == 0:
                return seg["vstart"]
            t = (pos - seg["gstart"]) / span
            return seg["vstart"] + t * (seg["vend"] - seg["vstart"])
    return segments[0]["vstart"] if pos < segments[0]["gstart"] else segments[-1]["vend"]


# ── Internal track builders ───────────────────────────────────────────────────

def _make_exon_track(
    exon_segs: list[dict],
    utr_segs: list[dict],
    segments: list[dict],
    total_vw: float,
    atg_pos: int,
    stop_pos: int,
    x_scale: alt.Scale,
    chart_width: int,
    exon_color: str,
    fontsize: int,
) -> alt.Chart:
    """Build the exon-structure layer (unconfigured, for composing)."""
    TRACK_H = 97
    MARKER_Y = 1      # top of ATG/Stop text
    ARROW_Y = 19      # tip of ATG/Stop triangle (pointing down toward exon)
    EXON_TOP = 32     # top of CDS rect
    EXON_BOT = 64     # bottom of CDS rect
    UTR_TOP = 37      # top of UTR rect
    UTR_BOT = 59      # bottom of UTR rect
    BACKBONE_Y = (EXON_TOP + EXON_BOT) // 2
    LABEL_Y = 76      # exon-number labels below

    _STROKE = "black"
    _SW = 0.8

    def _base(chart: alt.Chart) -> alt.Chart:
        return chart.properties(height=TRACK_H)

    layers: list[alt.Chart] = []

    # 1. Backbone line
    backbone_df = pd.DataFrame({"x": [0.0], "x2": [total_vw]})
    layers.append(_base(
        alt.Chart(backbone_df)
        .mark_rule(color="black", strokeWidth=1.5)
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            x2="x2:Q",
            y=alt.value(BACKBONE_Y),
        )
    ))

    # 2 & 3. Exon rectangles — sub-segments already split in _build_visual_segments
    cds_rows = [
        {"x": s["vstart"], "x2": s["vend"], "exon": f"Exon {s['exon']}"}
        for s in exon_segs
    ]
    utr_rows = [
        {"x": s["vstart"], "x2": s["vend"], "exon": f"UTR {s['exon']}"}
        for s in utr_segs
    ]

    if cds_rows:
        layers.append(_base(
            alt.Chart(pd.DataFrame(cds_rows))
            .mark_rect(fill=exon_color, stroke=_STROKE, strokeWidth=_SW)
            .encode(
                x=alt.X("x:Q", scale=x_scale, axis=None),
                x2="x2:Q",
                y=alt.value(EXON_TOP),
                y2=alt.value(EXON_BOT),
                tooltip=alt.Tooltip("exon:N", title="Region"),
            )
        ))
    if utr_rows:
        layers.append(_base(
            alt.Chart(pd.DataFrame(utr_rows))
            .mark_rect(fill=exon_color, stroke=_STROKE, strokeWidth=_SW)
            .encode(
                x=alt.X("x:Q", scale=x_scale, axis=None),
                x2="x2:Q",
                y=alt.value(UTR_TOP),
                y2=alt.value(UTR_BOT),
                tooltip=alt.Tooltip("exon:N", title="Region"),
            )
        ))

    # 4. Exon number labels (below) + "Exon" row label in left margin
    # Center each label over the full exon visual span (UTR + CDS sub-segments combined)
    exon_vspan: dict[str, list[float]] = {}
    for s in exon_segs + utr_segs:
        name = str(s["exon"])
        if name not in exon_vspan:
            exon_vspan[name] = [s["vstart"], s["vend"]]
        else:
            exon_vspan[name][0] = min(exon_vspan[name][0], s["vstart"])
            exon_vspan[name][1] = max(exon_vspan[name][1], s["vend"])
    label_df = pd.DataFrame([
        {
            "center": (span[0] + span[1]) / 2,
            # Strip a leading "X" from auto-generated names like "X1" → "1"
            "label": name[1:] if name[:1] == "X" and name[1:].isdigit() else name,
        }
        for name, span in sorted(exon_vspan.items(), key=lambda kv: kv[1][0])
    ])
    layers.append(_base(
        alt.Chart(label_df)
        .mark_text(fontSize=fontsize - 1, fontWeight="bold", baseline="top")
        .encode(
            x=alt.X("center:Q", scale=x_scale, axis=None),
            y=alt.value(LABEL_Y),
            text="label:N",
        )
    ))
    row_label_df = pd.DataFrame({"x": [-_LEFT_MARGIN / 2], "label": ["Exon"]})
    layers.append(_base(
        alt.Chart(row_label_df)
        .mark_text(fontSize=fontsize - 1, fontWeight="bold", baseline="top", align="center")
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            y=alt.value(LABEL_Y),
            text="label:N",
        )
    ))

    # 5. ATG / Stop: text label above, then triangle pointing down
    atg_vx = _gv(atg_pos, segments)
    stop_vx = _gv(stop_pos, segments)
    marker_df = pd.DataFrame({
        "x": [atg_vx, stop_vx],
        "label": ["ATG", "Stop"],
        "color": ["#1a7a1a", "#b22222"],
    })
    layers.append(_base(
        alt.Chart(marker_df)
        .mark_text(fontSize=fontsize - 2, fontWeight="bold", baseline="top", align="center")
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            y=alt.value(MARKER_Y),
            text="label:N",
            color=alt.Color("color:N", scale=None, legend=None),
        )
    ))
    layers.append(_base(
        alt.Chart(marker_df)
        .mark_point(shape="triangle-down", filled=True, size=fontsize * 7)
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            y=alt.value(ARROW_Y),
            color=alt.Color("color:N", scale=None, legend=None),
        )
    ))

    return alt.layer(*layers)


def _hex_darken(hex_color: str, factor: float) -> str:
    """Scale each RGB channel of a hex color by *factor* (0 = black, 1 = unchanged)."""
    c = hex_color.lstrip("#")
    r = min(255, int(int(c[0:2], 16) * factor))
    g = min(255, int(int(c[2:4], 16) * factor))
    b = min(255, int(int(c[4:6], 16) * factor))
    return f"#{r:02x}{g:02x}{b:02x}"


def _coverage_rects(lib_df: pd.DataFrame, segments: list[dict]) -> list[dict]:
    """Compute per-segment coverage depth in visual x-coordinates.

    Finds all visual x breakpoints (amplicon starts and ends), then counts
    how many amplicons overlap each consecutive pair.  Returns a list of
    dicts with keys ``x``, ``x2``, and ``depth`` (>= 1).
    """
    amplicons = [
        (_gv(int(row["start"]), segments), _gv(int(row["end"]), segments))
        for _, row in lib_df.iterrows()
    ]
    breakpoints = sorted({vx for vx_s, vx_e in amplicons for vx in (vx_s, vx_e)})

    rects: list[dict] = []
    for i in range(len(breakpoints) - 1):
        x0, x1 = breakpoints[i], breakpoints[i + 1]
        mid = (x0 + x1) / 2
        depth = sum(1 for vx_s, vx_e in amplicons if vx_s <= mid <= vx_e)
        if depth > 0:
            rects.append({"x": x0, "x2": x1, "depth": depth})
    return rects


def _covered_bases(lib_df: pd.DataFrame) -> int:
    """Total unique genomic bases covered by any library amplicon (union of intervals)."""
    intervals = sorted(zip(lib_df["start"].astype(int), lib_df["end"].astype(int)))
    merged: list[list[int]] = []
    for start, end in intervals:
        if merged and start < merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return sum(e - s for s, e in merged)


def _make_library_track(
    lib_df: pd.DataFrame,
    segments: list[dict],
    x_scale: alt.Scale,
    chart_width: int,
    lib_color: str,
    fontsize: int,
) -> alt.Chart:
    """Build the library-amplicon coverage track (unconfigured, for composing).

    All amplicons are drawn on a single horizontal band.  Where amplicons
    overlap the fill is progressively darkened: depth-1 regions use
    *lib_color* unchanged; each additional layer of coverage darkens the
    color by a fixed factor toward black.

    A bracket spanning the full covered region is drawn above the band.  The
    bracket opens in the middle to show the estimated variant count
    (``covered_bases × 3`` SNVs + ``covered_bases − 2`` 3-bp deletions)
    as a single "N Variants Designed" label.
    """
    # ── Layout ───────────────────────────────────────────────────────────
    BRACKET_Y = 20   # y of horizontal bracket arms
    V_PAD     = 32   # top of library band
    TICK_H    = V_PAD - BRACKET_Y - 5   # ticks stop 5px above the band top
    LANE_H    = 26
    BOT_PAD   = 8
    TRACK_H   = V_PAD + LANE_H + BOT_PAD
    _STROKE   = "black"
    _SW       = 0.8

    rects = _coverage_rects(lib_df, segments)
    if not rects:
        return alt.layer()

    def _base(chart: alt.Chart) -> alt.Chart:
        return chart.properties(height=TRACK_H)

    # ── Coverage rectangles ───────────────────────────────────────────────
    max_depth = max(r["depth"] for r in rects)

    def _color(depth: int) -> str:
        if max_depth == 1:
            return lib_color
        factor = 1.0 - (depth - 1) / max_depth * 0.65
        return _hex_darken(lib_color, factor)

    rect_df = pd.DataFrame([
        {"x": r["x"], "x2": r["x2"], "fill": _color(r["depth"]), "depth": r["depth"]}
        for r in rects
    ])
    layers: list[alt.Chart] = [
        _base(
            alt.Chart(rect_df)
            .mark_rect(stroke=_STROKE, strokeWidth=_SW)
            .encode(
                x=alt.X("x:Q", scale=x_scale, axis=None),
                x2="x2:Q",
                y=alt.value(V_PAD),
                y2=alt.value(V_PAD + LANE_H),
                color=alt.Color("fill:N", scale=None, legend=None),
                tooltip=alt.Tooltip("depth:Q", title="Coverage"),
            )
        )
    ]

    # ── Bracket above the library band ────────────────────────────────────
    cov = _covered_bases(lib_df)
    n_variants = cov * 3 + max(0, cov - 2)
    label = f"{n_variants:,} Variants Designed"

    all_vx = (
        [_gv(int(r["start"]), segments) for _, r in lib_df.iterrows()]
        + [_gv(int(r["end"]), segments) for _, r in lib_df.iterrows()]
    )
    x_lib_start = min(all_vx)
    x_lib_end   = max(all_vx)
    center      = (x_lib_start + x_lib_end) / 2
    gap_half    = len(label) * 4.0 + 10

    # Two horizontal arm segments (left and right of the text gap)
    arms_df = pd.DataFrame([
        {"x": x_lib_start, "x2": max(x_lib_start, center - gap_half)},
        {"x": min(x_lib_end, center + gap_half), "x2": x_lib_end},
    ])
    layers.append(_base(
        alt.Chart(arms_df)
        .mark_rule(color="black", strokeWidth=1)
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            x2="x2:Q",
            y=alt.value(BRACKET_Y),
        )
    ))

    # Vertical end-ticks at the outer ends of the bracket
    ticks_df = pd.DataFrame({"x": [x_lib_start, x_lib_end]})
    layers.append(_base(
        alt.Chart(ticks_df)
        .mark_rule(color="black", strokeWidth=1)
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            y=alt.value(BRACKET_Y),
            y2=alt.value(BRACKET_Y + TICK_H),
        )
    ))

    # Label centered in the gap, vertically centered on the bracket arms
    label_df = pd.DataFrame({"x": [center], "label": [label]})
    layers.append(_base(
        alt.Chart(label_df)
        .mark_text(fontSize=fontsize, fontWeight="bold", baseline="middle", align="center")
        .encode(
            x=alt.X("x:Q", scale=x_scale, axis=None),
            y=alt.value(BRACKET_Y),
            text="label:N",
        )
    ))

    # ── "N variant / libraries" row label in left margin ─────────────────
    n_libs = len(lib_df)
    mid_y = V_PAD + LANE_H // 2
    row_x = -_LEFT_MARGIN / 2
    for line, baseline in [(f"{n_libs} Variant", "bottom"), ("Libraries", "top")]:
        line_df = pd.DataFrame({"x": [row_x], "label": [line]})
        layers.append(_base(
            alt.Chart(line_df)
            .mark_text(fontSize=fontsize - 2, fontWeight="bold", baseline=baseline, align="center")
            .encode(
                x=alt.X("x:Q", scale=x_scale, axis=None),
                y=alt.value(mid_y),
                text="label:N",
            )
        ))

    return alt.layer(*layers)


# ── Public functions ──────────────────────────────────────────────────────────

def _parse_meta(metadata_df: pd.DataFrame, exon_df: pd.DataFrame) -> dict:
    meta = dict(zip(metadata_df["type"].str.lower(), metadata_df["info"]))
    strand = str(meta.get("strand", "plus")).strip().lower()
    # For minus-strand genes the 5′ end (ATG) is at the highest genomic
    # coordinate; adjust the fallback defaults accordingly.
    if strand == "minus":
        default_atg  = int(exon_df["end"].max())
        default_stop = int(exon_df["start"].min())
    else:
        default_atg  = int(exon_df["start"].min())
        default_stop = int(exon_df["end"].max())
    return {
        "atg_pos":    int(meta.get("atg", default_atg)),
        "stop_pos":   int(meta.get("stop", default_stop)),
        "exon_color": str(meta.get("exon_color", "#2E86C1")),
        "lib_color":  str(meta.get("lib_color", "#888888")),
        "strand":     strand,
    }


def _apply_strand(
    exon_df: pd.DataFrame,
    lib_df: pd.DataFrame | None,
    atg_pos: int,
    stop_pos: int,
    strand: str,
) -> tuple[pd.DataFrame, pd.DataFrame | None, int, int]:
    """Convert genomic coordinates to transcription-order space.

    For plus-strand genes all coordinates are returned unchanged.

    For minus-strand genes every coordinate is negated so that the 5′ end of
    the transcript (highest genomic position) becomes the most-negative value
    and therefore sorts leftmost.  Exon start/end are swapped after negation
    so that start < end is preserved in the transformed space.

    ``_build_visual_segments`` and ``_gv`` are strand-agnostic and work
    correctly on either original or negated coordinates.
    """
    if strand != "minus":
        return exon_df, lib_df, atg_pos, stop_pos

    exon_vgc = exon_df.copy()
    exon_vgc["start"] = -exon_df["end"].values
    exon_vgc["end"]   = -exon_df["start"].values

    lib_vgc: pd.DataFrame | None = None
    if lib_df is not None:
        lib_vgc = lib_df.copy()
        lib_vgc["start"] = -lib_df["end"].values
        lib_vgc["end"]   = -lib_df["start"].values

    return exon_vgc, lib_vgc, -atg_pos, -stop_pos


def make_exon_cartoon(
    exon_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    width: int = 800,
    intron_min_bp: int = 100,
    intron_max_bp: int = 10_000,
    intron_min_vw: float = 20.0,
    intron_max_vw: float = 60.0,
    utr_min_bp: int = 100,
    utr_max_bp: int = 5_000,
    utr_min_vw: float = 15.0,
    utr_max_vw: float = 50.0,
    fontsize: int = 16,
    exon_color: str | None = None,
) -> alt.Chart:
    """Draw a scalable exon-structure cartoon with compressed introns and UTRs.

    CDS exon widths are proportional to their genomic length.  Introns and
    UTR regions are each independently scaled between their respective
    *_min_vw* / *_max_vw* pixel limits based on genomic length, so very long
    UTRs or introns don't dominate the figure.

    A ``//`` marker appears at the midpoint of every intron.

    Colors are read from *metadata_df* (``exon_color`` row); the default is
    ``#2E86C1``.

    Args:
        exon_df: columns ``['exon', 'start', 'end']`` – genomic coordinates.
        metadata_df: columns ``['type', 'info']`` with rows ``'atg'``,
            ``'stop'``, and optionally ``'exon_color'``.
        width: visual width of the data area in pixels.
        intron_min_bp: introns shorter than this get *intron_min_vw* px.
        intron_max_bp: introns longer than this get *intron_max_vw* px.
        intron_min_vw: minimum intron visual width in pixels.
        intron_max_vw: maximum intron visual width in pixels.
        utr_min_bp: UTR regions shorter than this get *utr_min_vw* px.
        utr_max_bp: UTR regions longer than this get *utr_max_vw* px.
        utr_min_vw: minimum UTR visual width in pixels.
        utr_max_vw: maximum UTR visual width in pixels.
        fontsize: base font size for all text labels; individual elements
            scale relative to this value.
    """
    meta = _parse_meta(metadata_df, exon_df)
    if exon_color is not None:
        meta["exon_color"] = exon_color
    exon_df_vgc, _, atg_vgc, stop_vgc = _apply_strand(
        exon_df, None, meta["atg_pos"], meta["stop_pos"], meta["strand"]
    )
    segments, total_vw = _build_visual_segments(
        exon_df_vgc, width,
        intron_min_bp, intron_max_bp, intron_min_vw, intron_max_vw,
        utr_min_bp, utr_max_bp, utr_min_vw, utr_max_vw,
        atg_vgc, stop_vgc,
    )
    exon_segs = [s for s in segments if s["kind"] == "exon"]
    utr_segs = [s for s in segments if s["kind"] == "utr"]

    x_scale = alt.Scale(domain=[-_LEFT_MARGIN, total_vw])
    chart_width = width + _LEFT_MARGIN

    track = _make_exon_track(
        exon_segs, utr_segs, segments, total_vw,
        atg_vgc, stop_vgc,
        x_scale, chart_width, meta["exon_color"], fontsize,
    )
    return track.properties(width=chart_width).configure_axis(grid=False).configure_view(stroke=None)


def make_library_cartoon(
    exon_df: pd.DataFrame,
    lib_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    width: int = 800,
    intron_min_bp: int = 100,
    intron_max_bp: int = 10_000,
    intron_min_vw: float = 20.0,
    intron_max_vw: float = 60.0,
    utr_min_bp: int = 100,
    utr_max_bp: int = 5_000,
    utr_min_vw: float = 15.0,
    utr_max_vw: float = 50.0,
    fontsize: int = 16,
    exon_color: str | None = None,
    lib_color: str | None = None,
) -> alt.Chart:
    """Draw an exon-structure cartoon with a library-amplicon track below.

    The two tracks share the same compressed x-axis.  Overlapping amplicons
    are automatically stacked into lanes using a greedy left-to-right packing
    algorithm.

    Colors are read from *metadata_df* (``exon_color`` and ``lib_color``
    rows); defaults are ``#2E86C1`` and ``#888888`` respectively.

    Args:
        exon_df: columns ``['exon', 'start', 'end']`` – genomic coordinates.
        lib_df: columns ``['start', 'end']`` – amplicon genomic coordinates.
        metadata_df: columns ``['type', 'info']`` with rows ``'atg'``,
            ``'stop'``, ``'exon_color'``, and ``'lib_color'``.
        width: visual width of the data area in pixels.
        intron_min_bp / intron_max_bp / intron_min_vw / intron_max_vw:
            intron compression parameters (see ``make_exon_cartoon``).
        utr_min_bp / utr_max_bp / utr_min_vw / utr_max_vw:
            UTR compression parameters (see ``make_exon_cartoon``).
        fontsize: base font size for all text labels; individual elements
            scale relative to this value.
    """
    meta = _parse_meta(metadata_df, exon_df)
    if exon_color is not None:
        meta["exon_color"] = exon_color
    if lib_color is not None:
        meta["lib_color"] = lib_color
    exon_df_vgc, lib_df_vgc, atg_vgc, stop_vgc = _apply_strand(
        exon_df, lib_df, meta["atg_pos"], meta["stop_pos"], meta["strand"]
    )
    segments, total_vw = _build_visual_segments(
        exon_df_vgc, width,
        intron_min_bp, intron_max_bp, intron_min_vw, intron_max_vw,
        utr_min_bp, utr_max_bp, utr_min_vw, utr_max_vw,
        atg_vgc, stop_vgc,
    )
    exon_segs = [s for s in segments if s["kind"] == "exon"]
    utr_segs = [s for s in segments if s["kind"] == "utr"]

    x_scale = alt.Scale(domain=[-_LEFT_MARGIN, total_vw])
    chart_width = width + _LEFT_MARGIN

    exon_track = _make_exon_track(
        exon_segs, utr_segs, segments, total_vw,
        atg_vgc, stop_vgc,
        x_scale, chart_width, meta["exon_color"], fontsize,
    )
    lib_track = _make_library_track(
        lib_df_vgc, segments, x_scale, chart_width, meta["lib_color"], fontsize,
    )

    return (
        alt.vconcat(exon_track, lib_track, spacing=8)
        .properties(width=chart_width)
        .configure_axis(grid=False)
        .configure_view(stroke=None)
    )
