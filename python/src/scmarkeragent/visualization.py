#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render reference-free production figures from language-neutral figure data.

The shared reporter writes the data once, so Python and R annotation arms render
the same production outputs without re-annotation:

  1 umap_cluster      UMAP colored by cluster id, labels on data
  2 umap_celltype     UMAP colored by primary_annotation (finest call), legend = "cluster: primary"
  3 dotplot_celltype_markers  LLM-validated, source-backed identity markers
                              y-axis = full primary labels; x-block labels may truncate

Style = the user's established EvidenCellMarker look (eviden 25-color palette;
dotplot as shape-21 black-stroke points sized by pct_exp, filled by avg_exp_scaled
on the red gradient #fee0d2 -> #fc9272 -> #de2d26) + academic-figure-skill export
rigor (Arial/Liberation Sans, clean spines, vector PDF + 300 dpi PNG, print dims).

This is a library: `plot_dataset` is the entry point, and nothing here is a script.
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from . import plot_params as pp  # noqa: E402

FIGDATA_DIR = "figure_data"
FIGURES_DIR = "figures"

# ---- eviden palette defaults (verbatim from EvidenCellMarker/.../visualization.R) -----
_EVIDEN_DEFAULT = [
    "#a4cde1",
    "#67a4cc",
    "#277fb8",
    "#549da3",
    "#96cb8f",
    "#8bc96d",
    "#4dae47",
    "#5c9e43",
    "#b79973",
    "#f38989",
    "#ec5051",
    "#e32427",
    "#ef6a45",
    "#f9b769",
    "#f9a341",
    "#f48521",
    "#ee8e46",
    "#d4a6a8",
    "#af93c4",
    "#8660a8",
    "#815e99",
    "#c6b598",
    "#f6f28f",
    "#d4a55b",
    "#b05a28",
]

# ---- resolve the shared figure params (EDIT config/plot_params.json, NOT this file) ---
_FONT_FAMILY = pp.text("global", "font_family", "Arial")
_FONT_FALLBACK = pp.str_list(
    "global", "font_fallback", ["Helvetica", "Liberation Sans", "DejaVu Sans"]
)
_BASE_FS = pp.num("global", "base_font_size", 8, lo=1)
_DPI = pp.num("global", "dpi", 300, lo=1)
_FORMATS = [f.lower() for f in pp.str_list("global", "formats", ["pdf", "png"])]
_FIG_BG = pp.text("global", "figure_background", "white")

# ---- academic-figure-skill typography + export baseline (Arial -> Liberation Sans) ----
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": [_FONT_FAMILY] + _FONT_FALLBACK,
        "font.size": _BASE_FS,
        "axes.titlesize": _BASE_FS,
        "axes.labelsize": _BASE_FS,
        "xtick.labelsize": _BASE_FS - 1,
        "ytick.labelsize": _BASE_FS - 1,
        "legend.fontsize": _BASE_FS - 1,
        "figure.titlesize": _BASE_FS + 1,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": _DPI,
        "figure.facecolor": _FIG_BG,
        "savefig.facecolor": _FIG_BG,
    }
)
MM = 1 / 25.4  # inches per millimetre (matplotlib figure/inset sizes are in inches)

EVIDEN_BASE = pp.str_list("palette", "eviden_base", _EVIDEN_DEFAULT)
# eviden dotplot red gradient (fill mapped on avg.exp.scaled), 3-stop low->mid->high
DOT_CMAP = LinearSegmentedColormap.from_list(
    "eviden_red",
    pp.str_list("dotplot", "gradient_colors", ["#fee0d2", "#fc9272", "#de2d26"]),
)


def eviden_palette(n):
    """n distinct colors: eviden base first, then tab20/tab20b hues, then HSV fallback
    (mirrors .eviden_build_palette semantics: user colors, then extend, then distinct)."""
    if n <= len(EVIDEN_BASE):
        return EVIDEN_BASE[:n]
    pal = list(EVIDEN_BASE)
    for cmap in ("tab20", "tab20b", "tab20c"):
        cols = [mpl.colors.to_hex(c) for c in plt.get_cmap(cmap).colors]
        pal += [c for c in cols if c not in pal]
        if len(pal) >= n:
            break
    if len(pal) < n:
        extra = plt.get_cmap("hsv")(np.linspace(0, 1, n - len(pal), endpoint=False))
        pal += [mpl.colors.to_hex(c) for c in extra]
    return pal[:n]


def _save(fig, stem):
    for fmt in _FORMATS:
        fig.savefig(f"{stem}.{fmt}", bbox_inches="tight", dpi=_DPI, facecolor=_FIG_BG)
    plt.close(fig)


def _point_size(n):
    fixed = pp.num("umap", "point_size", 0, lo=0)
    if fixed and fixed > 0:
        return float(fixed)
    return max(1.5, min(14.0, 90000.0 / max(n, 1)))


def _area_to_diam(area):
    """matplotlib scatter s is AREA (pt^2); Line2D markersize is DIAMETER (pt)."""
    return 2.0 * np.sqrt(max(float(area), 0.0) / np.pi)


def _dot_area(pct, smin, smax, cap):
    """Dot AREA (pt^2) for a pct_exp, following ggplot2 area_pal (the Seurat/eviden
    convention rflow/plots.R uses): rendered area =
        ( sqrt(smin) + (sqrt(smax)-sqrt(smin)) * sqrt(pct/100) )^2 , clamped to cap.
    Endpoints give exactly smin (pct=0) and smax (pct=100); the sqrt interior makes the
    Python dots match R's scale_size() dots at EVERY pct, not just the endpoints."""
    v = np.clip(np.asarray(pct, dtype=float) / 100.0, 0.0, 1.0)
    a = (np.sqrt(smin) + (np.sqrt(smax) - np.sqrt(smin)) * np.sqrt(v)) ** 2
    return np.minimum(a, cap)


# ------------------------------------------------------------------ UMAP panels ----
# Academic-figure consistency rule (user requirement): every UMAP's DATA panel must be an
# identical fixed-size SQUARE with the SAME shared data limits + equal aspect, so the
# umap_cluster / umap_celltype / umap_reference panels (and the same panel across datasets)
# are directly comparable side by side. A category legend of variable size would otherwise
# distort the data area, so the legend is rendered as a SEPARATE panel and then STITCHED to
# the right of the fixed data panel (never resizing it). We therefore never use tight bbox
# for UMAPs: the figure canvas is a fixed physical size with fixed margins.


def _save_fixed(fig, stem):
    """Save a UMAP figure at its EXACT set size (no tight bbox), so the data axes keeps an
    identical physical size across every panel and dataset (comparability guarantee)."""
    for fmt in _FORMATS:
        fig.savefig(f"{stem}.{fmt}", dpi=_DPI, facecolor=_FIG_BG)
    plt.close(fig)


def umap_limits(cells):
    """Shared SQUARE data limits for ALL UMAP panels of a dataset: one common bounding box
    (centered, equal x/y span with a small margin) computed once from every cell, so panels
    that show only a subset (e.g. reference-labelled cells) sit in the identical coordinate
    frame as the full-cell panels."""
    x, y = cells["umap_x"].astype(float), cells["umap_y"].astype(float)
    cx, cy = (x.min() + x.max()) / 2.0, (y.min() + y.max()) / 2.0
    half = 0.5 * max(x.max() - x.min(), y.max() - y.min()) * 1.05 or 1.0
    return (cx - half, cx + half), (cy - half, cy + half)


def _umap_panel(
    cells,
    cats,
    color_of,
    title,
    stem,
    on_data_labels=False,
    order=None,
    xlim=None,
    ylim=None,
):
    """Render ONE fixed-size square UMAP data panel (no legend). Same physical size and, when
    xlim/ylim are passed, the same data frame for every panel -> directly comparable."""
    n = len(cells)
    alpha = pp.num("umap", "point_alpha", 1.0, lo=0, hi=1)
    rasterize = pp.boolean("umap", "rasterize", True)
    label_fs = pp.num("umap", "label_font_size", 7, lo=1)
    panel_mm = pp.num("umap", "panel_mm", 55, lo=10)  # square DATA-area side (mm)
    m_l = pp.num("umap", "margin_left_mm", 8, lo=0)
    m_r = pp.num("umap", "margin_right_mm", 3, lo=0)
    m_b = pp.num("umap", "margin_bottom_mm", 7, lo=0)
    m_t = pp.num("umap", "margin_top_mm", 9, lo=0)
    W = (panel_mm + m_l + m_r) * MM
    H = (panel_mm + m_b + m_t) * MM
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes(
        [(m_l * MM) / W, (m_b * MM) / H, (panel_mm * MM) / W, (panel_mm * MM) / H]
    )
    ps = _point_size(n)
    cat_list = order if order is not None else sorted(pd.unique(cats))
    for cat in cat_list:
        m = cats == cat
        if not m.any():
            continue
        ax.scatter(
            cells.loc[m, "umap_x"],
            cells.loc[m, "umap_y"],
            s=ps,
            c=color_of[cat],
            linewidths=0,
            alpha=alpha,
            rasterized=rasterize,
            label=str(cat),
        )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(title, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_aspect("equal")  # no distortion; square framing
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.6)
    if on_data_labels:
        for cat in cat_list:
            m = cats == cat
            if m.any():
                ax.text(
                    cells.loc[m, "umap_x"].median(),
                    cells.loc[m, "umap_y"].median(),
                    str(cat),
                    fontsize=label_fs,
                    fontweight="bold",
                    ha="center",
                    va="center",
                    color="black",
                    bbox=dict(
                        boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.65
                    ),
                )
    _save_fixed(fig, stem)


def _umap_legend(order, color_of, legend_title, stem):
    """Render the category legend as its OWN standalone panel (keys = filled dots), sized to the
    entries. Kept separate from the data panel so the data panel size never changes; the two are
    stitched by _stitch_umap for the delivered '<stem>_labeled.png'."""
    leg_ms = pp.num("umap", "legend_marker_size", 4, lo=0.1)
    leg_ncol = int(pp.num("umap", "legend_ncol", 1, lo=1))
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=leg_ms,
            markerfacecolor=color_of[c],
            markeredgewidth=0,
            label=str(c),
        )
        for c in order
    ]
    n_rows = int(np.ceil(len(order) / max(leg_ncol, 1)))
    longest = max((len(str(c)) for c in order), default=8)
    w_in = (10 + leg_ncol * (6 + longest * 1.5)) * MM
    h_in = max(20, 6 + n_rows * 4.4) * MM
    fig = plt.figure(figsize=(w_in, h_in))
    fig.legend(
        handles=handles,
        title=legend_title,
        loc="center left",
        frameon=False,
        ncol=leg_ncol,
        fontsize=_BASE_FS - 2,
        title_fontsize=_BASE_FS - 1,
        handletextpad=0.3,
        labelspacing=0.35,
        borderaxespad=0.2,
    )
    for fmt in _FORMATS:
        fig.savefig(f"{stem}.{fmt}", dpi=_DPI, facecolor=_FIG_BG)
    plt.close(fig)


def _stitch_umap(panel_png, legend_png, out_png):
    """Horizontally stitch the fixed data panel + its legend into one PNG (data panel left,
    legend right, vertically centered). The data panel keeps its exact size, so the stitched
    figure stays comparable across panels/datasets."""
    try:
        from PIL import Image
    except Exception:
        return
    if not (os.path.exists(panel_png) and os.path.exists(legend_png)):
        return
    a, b = Image.open(panel_png).convert("RGBA"), Image.open(legend_png).convert("RGBA")
    h = max(a.height, b.height)
    canvas = Image.new("RGBA", (a.width + b.width, h), (255, 255, 255, 255))
    canvas.paste(a, (0, (h - a.height) // 2), a)
    canvas.paste(b, (a.width, (h - b.height) // 2), b)
    canvas.convert("RGB").save(out_png, dpi=(_DPI, _DPI))


def _umap(
    cells,
    cats,
    color_of,
    title,
    stem,
    on_data_labels=False,
    legend_title=None,
    order=None,
    xlim=None,
    ylim=None,
):
    """Render a UMAP: always the fixed square data panel (<stem>); when it carries a category
    legend, also render the standalone legend (<stem>_legend) and the stitched combined PNG
    (<stem>_labeled.png) so the deliverable is 'drawn separately then stitched' per requirement."""
    cat_list = order if order is not None else sorted(pd.unique(cats))
    _umap_panel(
        cells,
        cats,
        color_of,
        title,
        stem,
        on_data_labels=on_data_labels,
        order=cat_list,
        xlim=xlim,
        ylim=ylim,
    )
    if (not on_data_labels) and pp.boolean("umap", "separate_legend", True):
        _umap_legend(cat_list, color_of, legend_title, stem + "_legend")
        _stitch_umap(stem + ".png", stem + "_legend.png", stem + "_labeled.png")


def _axes_size_mm(fig, ax):
    """The data area's own size in mm. Anything placed beside the panel is measured against
    this and never against the figure, whose margins are a fraction of a size that itself
    varies with the number of genes and rows."""
    figure_width, figure_height = fig.get_size_inches()
    box = ax.get_position()
    return box.width * figure_width / MM, box.height * figure_height / MM


# colorbar side placement: (loc, x-anchor for vertical bar) relative to the axes
def _ensure_legend_gutter(fig, ax, position, reserve_mm):
    """Ensure ≥`reserve_mm` of figure space is free on `position` of `ax`.

    Shrinks the data axes only when the current margin cannot hold the legend column
    (narrow panels). Wide panels already have spare margin; `_legend_anchor` then
    parks the legend in that margin instead of hugging the data at axes+epsilon."""
    if reserve_mm <= 0:
        return
    fw_mm = fig.get_size_inches()[0] / MM
    fh_mm = fig.get_size_inches()[1] / MM
    pos = ax.get_position()
    min_span = 0.35
    if position == "right":
        need = reserve_mm / fw_mm
        if 1.0 - pos.x1 + 1e-9 >= need:
            return
        new_x1 = 1.0 - need
        if new_x1 - pos.x0 >= min_span:
            ax.set_position([pos.x0, pos.y0, new_x1 - pos.x0, pos.height])
    elif position == "left":
        need = reserve_mm / fw_mm
        if pos.x0 + 1e-9 >= need:
            return
        if pos.x1 - need >= min_span:
            ax.set_position([need, pos.y0, pos.x1 - need, pos.height])
    elif position == "top":
        need = reserve_mm / fh_mm
        if 1.0 - pos.y1 + 1e-9 >= need:
            return
        new_y1 = 1.0 - need
        if new_y1 - pos.y0 >= min_span:
            ax.set_position([pos.x0, pos.y0, pos.width, new_y1 - pos.y0])
    elif position == "bottom":
        need = reserve_mm / fh_mm
        if pos.y0 + 1e-9 >= need:
            return
        if pos.y1 - need >= min_span:
            ax.set_position([pos.x0, need, pos.width, pos.y1 - need])


def _legend_anchor(fig, ax, position, gap_mm, column_mm, max_gap_mm=20.0):
    """Axes-fraction anchor for the legend column's near edge.

    Failure mode this fixes: legend at axes+epsilon while unused figure margin (and
    `bbox_inches=tight` overhang from rotated block labels) sits *beyond* the legend,
    so the colorbar looks glued onto the last gene columns — especially on wide, short
    by-celltype panels. Park the column in that gutter, but never more than
    `max_gap_mm` from the data so wide panels do not grow a huge empty band."""
    fw_mm = fig.get_size_inches()[0] / MM
    fh_mm = fig.get_size_inches()[1] / MM
    pos = ax.get_position()
    width_mm = max(pos.width * fw_mm, 1e-9)
    height_mm = max(pos.height * fh_mm, 1e-9)
    if position == "right":
        gutter = (1.0 - column_mm / fw_mm - pos.x0) / max(pos.width, 1e-9)
        lo = 1.0 + gap_mm / width_mm
        hi = 1.0 + max_gap_mm / width_mm
        return min(max(gutter, lo), hi), 1.0
    if position == "left":
        gutter = (column_mm / fw_mm - pos.x0) / max(pos.width, 1e-9)
        lo = -max_gap_mm / width_mm
        hi = -gap_mm / width_mm
        return max(min(gutter, hi), lo), 1.0
    if position == "top":
        gutter = (1.0 - column_mm / fh_mm - pos.y0) / max(pos.height, 1e-9)
        lo = 1.0 + gap_mm / height_mm
        hi = 1.0 + max_gap_mm / height_mm
        return 0.0, min(max(gutter, lo), hi)
    if position == "bottom":
        gutter = (column_mm / fh_mm - pos.y0) / max(pos.height, 1e-9)
        lo = -max_gap_mm / height_mm
        hi = -gap_mm / height_mm
        return 0.0, max(min(gutter, hi), lo)
    return 1.0 + gap_mm / width_mm, 1.0


def _colorbar_inset(fig, sc, ax, orient, position, length_mm, thick_mm, gap_mm, column_mm):
    """Absolute-size (mm) colorbar via native `Axes.inset_axes`.

    Must NOT use `mpl_toolkits.axes_grid1.inset_locator.inset_axes`: with
    `bbox_inches='tight'`, the PDF backend displaces the colorbar's rasterized
    gradient mesh while leaving the frame/ticks in place (matplotlib#27143 /
    #27763) — PNG is fine, which is exactly the PDF-only failure we hit. Native
    `ax.inset_axes` keeps mesh and frame aligned under tight cropping.

    Anchored via `_legend_anchor` so the bar sits in the figure gutter rather than
    flush on the data. Returns (colorbar, occupied_fraction) for placing the size
    legend clear of the bar."""
    width_mm, height_mm = _axes_size_mm(fig, ax)
    if orient == "horizontal":
        w_frac, h_frac = length_mm / width_mm, thick_mm / height_mm
    else:
        w_frac, h_frac = thick_mm / width_mm, length_mm / height_mm
    ax_x, ax_y = _legend_anchor(fig, ax, position, gap_mm, column_mm)
    # bounds = [x0, y0, width, height] in axes fraction; near-edge from _legend_anchor.
    if position == "right":
        bounds = [ax_x, 1.0 - h_frac, w_frac, h_frac]
    elif position == "left":
        bounds = [ax_x - w_frac, 1.0 - h_frac, w_frac, h_frac]
    elif position == "top":
        bounds = [0.0, ax_y, w_frac, h_frac]
    else:  # bottom
        bounds = [0.0, ax_y - h_frac, w_frac, h_frac]
    cax = ax.inset_axes(bounds)
    cb = fig.colorbar(sc, cax=cax, orientation=orient)
    # Keep ticks on the outer side of the bar so numbers sit in the gutter, not
    # in the gap between data and colorbar (matplotlib sometimes flips them).
    if orient == "vertical" and position == "right":
        cb.ax.yaxis.set_ticks_position("right")
        cb.ax.yaxis.set_label_position("right")
    elif orient == "vertical" and position == "left":
        cb.ax.yaxis.set_ticks_position("left")
        cb.ax.yaxis.set_label_position("left")
    return cb, length_mm / (height_mm if orient == "vertical" else width_mm)


# --------------------------------------------------------------- dotplot panel ----
def _dot_gradient_mid():
    """Representative fill for the pct.exp size-legend keys = midpoint of the fill gradient
    (so each key is ONE shape-21 dot: gradient fill + concentric black stroke, not a hollow
    ring nor an offset fill/border pair)."""
    override = pp.text("dotplot", "legend_dot_fill", "")
    if override:
        return override
    return mpl.colors.to_hex(DOT_CMAP(0.5))


def _dotplot(
    long,
    title,
    stem,
    grouped=False,
    y_label="Cluster: cell type",
    row_column="cluster_celltype",
    row_order_column="cluster_order",
    size_column="pct_exp",
    size_label="pct.exp",
):
    """Unified dotplot renderer. When `grouped` and the long matrix carries gene_group /
    gene_group_order, genes are laid out BLOCK-DIAGONAL: each cell type's markers form a
    contiguous, separated, labeled block (the EvidenCellMarker plot_marker_dotplot look).
    The size legend is drawn as single shape-21 dots (gradient fill + black stroke).

    Rows and the size channel are parameters because the same panel is delivered per
    cluster and per resolved cell type, and in a publication-support variant where the dot
    area is the marker's curated publication count instead of its detection percentage.
    Keeping one renderer is what stops those variants from drifting apart visually."""
    if long is None or not len(long):
        print(f"  [plots] dotplot '{title}': no data; skip")
        return
    smin, smax = pp.num_list("dotplot", "dot_size_range", [1.0, 55.0], n=2, lo=0)
    cap = pp.num("dotplot", "dot_size_max_cap", 55.0, lo=0)
    stroke = pp.num("dotplot", "dot_stroke", 0.9, lo=0)
    cb_len = pp.num("dotplot", "colorbar_length_mm", 18, lo=1)
    cb_thick = pp.num("dotplot", "colorbar_thickness_mm", 3.0, lo=0.1)
    cb_orient = pp.text(
        "dotplot",
        "colorbar_orientation",
        "vertical",
        choices={"vertical", "horizontal"},
    )
    cb_pos = pp.text(
        "dotplot",
        "colorbar_position",
        "right",
        choices={"right", "left", "top", "bottom"},
    )
    legend_gap = pp.num("dotplot", "legend_gap_mm", 4.0, lo=0)
    tick_fs = pp.num("dotplot", "tick_font_size", 6, lo=1)
    x_angle = pp.num("dotplot", "x_axis_angle", 90)
    grid_on = pp.boolean("dotplot", "panel_grid", True)
    grid_alpha = pp.num("dotplot", "panel_grid_alpha", 0.3, lo=0, hi=1)
    max_genes = int(pp.num("dotplot", "max_genes", 60, lo=0))
    w_per = pp.num("dotplot", "width_mm_per_gene", 3.4)
    w_base = pp.num("dotplot", "width_mm_base", 55)
    w_min = pp.num("dotplot", "width_mm_min", 120)
    h_per = pp.num("dotplot", "height_mm_per_cluster", 6.0)
    h_base = pp.num("dotplot", "height_mm_base", 32)
    h_min = pp.num("dotplot", "height_mm_min", 60)
    # block-diagonal grouping knobs (evidence panels)
    has_group = (
        grouped
        and ("gene_group" in long.columns)
        and (long["gene_group"].astype(str) != "").any()
    )
    block_gap = pp.num("dotplot", "block_gap", 0.9, lo=0)
    sep_on = pp.boolean("dotplot", "block_separator", True)
    sep_alpha = pp.num("dotplot", "block_separator_alpha", 0.55, lo=0, hi=1)
    lab_on = pp.boolean("dotplot", "show_block_labels", True)
    lab_fs = pp.num("dotplot", "block_label_font_size", 5.5, lo=1)
    lab_ang = pp.num("dotplot", "block_label_angle", 30)
    lab_max = int(pp.num("dotplot", "block_label_max_chars", 22, lo=1))

    slot_column = "marker_slot" if "marker_slot" in long.columns else "gene"
    slot_rows = long.sort_values("gene_order").drop_duplicates(slot_column)
    slots = slot_rows[slot_column].tolist()
    slot_labels = dict(zip(slot_rows[slot_column], slot_rows["gene"]))
    if max_genes and len(slots) > max_genes:
        slots = slots[:max_genes]
        long = long[long[slot_column].isin(set(slots))].copy()
    grp = long.sort_values(row_order_column)[row_column].drop_duplicates().tolist()
    yi = {c: i for i, c in enumerate(grp)}
    # The size channel is scaled against its own maximum so a count-valued channel and a
    # percentage one share the area_pal curve. A percentage is spread evenly over its range
    # and is scaled directly; publication counts are not -- they run from one paper to
    # several hundred, so a linear scale would collapse every ordinary marker to the
    # minimum dot and show only the handful of famous genes. Counts are therefore placed on
    # a log scale, which is also the scale on which one more paper is worth reading as a
    # difference. The legend keys are drawn from the same scale, so a reader reads the dots
    # off the counts they were built from.
    size_values = np.asarray(long[size_column].values, dtype=float)
    size_max = float(np.nanmax(size_values)) if size_values.size else 0.0
    if not np.isfinite(size_max) or size_max <= 0:
        size_max = 1.0
    if size_column == "pct_exp":

        def _size_fraction(value):
            return np.asarray(value, dtype=float) / size_max

        size_keys = [0, 25, 50, 75, 100]
    else:

        def _size_fraction(value):
            counts = np.clip(np.asarray(value, dtype=float), 0.0, None)
            return np.log1p(counts) / np.log1p(size_max)

        size_keys = sorted(
            {
                int(round(np.expm1(f * np.log1p(size_max))))
                for f in (0.0, 0.25, 0.5, 0.75, 1.0)
            }
        )
    size_fraction = _size_fraction(size_values)
    # per-gene x coordinate with a gap between blocks (block index from gene_group_order)
    gblock, glabel = {}, {}
    if has_group:
        gg = long.drop_duplicates(slot_column)[
            [slot_column, "gene_group", "gene_group_order"]
        ]
        gblock = dict(zip(gg[slot_column], gg["gene_group_order"].astype(int)))
        glabel = dict(zip(gg[slot_column], gg["gene_group"].astype(str)))
    xpos = {}
    for i, slot in enumerate(slots):
        xpos[slot] = i + (block_gap * int(gblock.get(slot, 0)) if has_group else 0.0)
    xmax = max(xpos.values()) if xpos else 0.0
    vmin, vmax = long["avg_exp_scaled"].min(), long["avg_exp_scaled"].max()
    if vmin == vmax:
        vmin, vmax = vmin - 1e-6, vmax + 1e-6

    h_mm = max(h_min, len(grp) * h_per + h_base)
    w = max(w_min, (xmax + 1) * w_per + w_base) * MM
    h = h_mm * MM
    fig, ax = plt.subplots(figsize=(w, h))
    x = long[slot_column].map(xpos).values
    y = long[row_column].map(yi).values
    sizes = _dot_area(100.0 * size_fraction, smin, smax, cap)  # area_pal curve
    sc = ax.scatter(
        x,
        y,
        s=sizes,
        c=long["avg_exp_scaled"].values,
        cmap=DOT_CMAP,
        vmin=vmin,
        vmax=vmax,
        marker="o",
        edgecolors="black",
        linewidths=stroke,
    )
    ax.set_xticks([xpos[slot] for slot in slots])
    ax.set_xticklabels(
        [slot_labels[slot] for slot in slots],
        rotation=x_angle,
        fontsize=tick_fs,
        ha="center",
        va="top",
    )
    ax.set_yticks(range(len(grp)))
    # Y-axis shows the full primary annotation (no truncation). Block labels drawn
    # above the x-axis may still be shortened via block_label_max_chars.
    ax.set_yticklabels(grp, fontsize=tick_fs)
    ax.set_xlim(-0.6, xmax + 0.6)
    # Keep the visual row order identical to cluster_summary.csv: cluster 0 at
    # the top, then ascending cluster order downward.
    ax.set_ylim(len(grp) - 0.4, -0.6)
    ax.set_xlabel("Markers")
    ax.set_ylabel(y_label)
    # When per-block cell-type labels are drawn above the axes (grouped evidence panels),
    # push the title up by their worst-case vertical rise so it never overlaps them.
    title_pad = mpl.rcParams.get("axes.titlepad", 6.0)
    if has_group and lab_on:
        title_pad += abs(np.sin(np.radians(lab_ang))) * lab_max * lab_fs * 0.6 + lab_fs
    ax.set_title(title, fontweight="bold", pad=title_pad)
    if grid_on:
        ax.grid(True, color="grey", linewidth=0.3, alpha=grid_alpha)
    ax.set_axisbelow(True)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.6)
    # block-diagonal separators + labels (evidence panels only)
    if has_group:
        blocks = {}
        for slot in slots:
            blocks.setdefault(gblock.get(slot, 0), []).append(xpos[slot])
        ordered = sorted(blocks)
        for bi in ordered[:-1]:
            xr = max(blocks[bi])
            xnext = min(blocks[ordered[ordered.index(bi) + 1]])
            if sep_on:
                ax.axvline(
                    (xr + xnext) / 2.0, color="grey", linewidth=0.4, alpha=sep_alpha
                )
        if lab_on:
            for bi in ordered:
                xs = blocks[bi]
                lab = next(
                    (glabel[slot] for slot in slots if gblock.get(slot, 0) == bi),
                    "",
                )
                lab = lab if len(lab) <= lab_max else lab[: lab_max - 1] + "\u2026"
                ax.text(
                    (min(xs) + max(xs)) / 2.0,
                    -0.3,
                    lab,
                    rotation=lab_ang,
                    fontsize=lab_fs,
                    ha="left" if lab_ang else "center",
                    va="bottom",
                    rotation_mode="anchor",
                    color="black",
                )
    # Legend column: clearance + bar + tick/label/size-key width. Reserve gutter on
    # narrow panels; on wide panels park the column in the existing right margin
    # (capped) so it does not sit flush on the last gene columns.
    legend_col_mm = legend_gap + cb_thick + 14.0
    _ensure_legend_gutter(fig, ax, cb_pos, legend_col_mm)
    cb, cb_fraction = _colorbar_inset(
        fig, sc, ax, cb_orient, cb_pos, cb_len, cb_thick, legend_gap, legend_col_mm
    )
    cb.set_label("avg.exp.scaled", fontsize=_BASE_FS - 1)
    cb.ax.tick_params(labelsize=tick_fs)
    cb.outline.set_linewidth(0.4)
    # pct.exp size legend: single shape-21 dots (gradient fill + concentric black stroke),
    # placed just below the (short) colorbar on the same side. FIX: one Line2D per key with
    # markerfacecolor (gradient) AND markeredgecolor="black" on the SAME marker -- so the fill
    # and the black border are one aligned dot, never a hollow ring or an offset fill/border.
    # The colorbar's own share of the axes is what the legend clears, so a panel with few
    # rows cannot place the keys on top of the bar.
    leg_fill = _dot_gradient_mid()
    axes_width_mm, axes_height_mm = _axes_size_mm(fig, ax)
    leg_x, _ = _legend_anchor(fig, ax, cb_pos, legend_gap, legend_col_mm)
    leg_y = (
        (1.0 - cb_fraction - legend_gap / axes_height_mm)
        if cb_orient == "vertical" and cb_pos in ("right", "left")
        else 0.5
    )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=leg_fill,
            markeredgecolor="black",
            markeredgewidth=stroke,
            markersize=_area_to_diam(
                _dot_area(100.0 * _size_fraction(key), smin, smax, cap)
            ),
            label=f"{key}",
        )
        for key in size_keys
    ]
    ax.legend(
        handles=handles,
        title=size_label,
        loc="upper left",
        bbox_to_anchor=(leg_x, leg_y),
        fontsize=tick_fs,
        title_fontsize=_BASE_FS - 1,
        labelspacing=0.7,
        borderaxespad=0,
    )
    _save(fig, stem)


# ------------------------------------------------------------------- driver ----
def plot_dataset(tag, outdir, fig_subdir=FIGURES_DIR):
    """Render reference-free figures for one completed production run."""
    os.makedirs(outdir, exist_ok=True)
    bundle = os.path.join(outdir, FIGDATA_DIR)
    figdir = os.path.join(outdir, fig_subdir) if fig_subdir else outdir
    os.makedirs(figdir, exist_ok=True)

    cells = pd.read_csv(
        os.path.join(bundle, "cells.csv"), dtype={"cluster": str}, keep_default_na=False
    )
    cm = pd.read_csv(
        os.path.join(bundle, "clustermap.csv"),
        dtype={"cluster": str},
        keep_default_na=False,
    )
    identity_markers = pd.read_csv(os.path.join(bundle, "dotplot_celltype_markers.csv"))
    celltype_marker_path = os.path.join(
        bundle, "dotplot_celltype_markers_by_celltype.csv"
    )
    celltype_markers = (
        pd.read_csv(celltype_marker_path)
        if os.path.isfile(celltype_marker_path)
        else None
    )
    if "umap_x" not in cells or cells["umap_x"].isna().all():
        print(f"  [plots] {tag}: no UMAP coords in bundle; skipping UMAP panels")
        umap_ok = False
    else:
        umap_ok = True

    order = cm.sort_values("cluster_order")
    clusters = order["cluster"].tolist()
    cc_of = dict(zip(order["cluster"], order["cluster_celltype"]))
    pal = eviden_palette(len(clusters))
    color_by_cluster = {c: pal[i] for i, c in enumerate(clusters)}
    color_by_cc = {cc_of[c]: pal[i] for i, c in enumerate(clusters)}
    cc_order = [cc_of[c] for c in clusters]

    if umap_ok:
        # ONE shared square data frame for every UMAP panel of this dataset (comparability).
        xlim, ylim = umap_limits(cells)
        _umap(
            cells,
            cells["cluster"],
            color_by_cluster,
            "Clusters",
            os.path.join(figdir, "umap_cluster"),
            on_data_labels=True,
            order=clusters,
            xlim=xlim,
            ylim=ylim,
        )
        cells_cc = cells.assign(cc=cells["cluster"].map(cc_of))
        _umap(
            cells_cc,
            cells_cc["cc"],
            color_by_cc,
            "Primary annotation",
            os.path.join(figdir, "umap_celltype"),
            on_data_labels=False,
            legend_title="Cluster: primary annotation",
            order=cc_order,
            xlim=xlim,
            ylim=ylim,
        )
    _dotplot(
        identity_markers,
        "LLM-validated cell-type identity markers",
        os.path.join(figdir, "dotplot_celltype_markers"),
        grouped=True,
        y_label="Cluster: primary annotation",
    )
    if celltype_markers is not None:
        # The same markers with clusters sharing a label pooled into one row, plus a
        # support-weighted twin whose dot area is the curated publication count, so a
        # reader can see identity evidence and literature backing side by side.
        _dotplot(
            celltype_markers,
            "Cell-type identity markers",
            os.path.join(figdir, "dotplot_celltype_markers_by_celltype"),
            grouped=True,
            y_label="Cell type",
            row_column="cell_type",
            row_order_column="cell_type_order",
        )
        _dotplot(
            celltype_markers,
            "Cell-type identity markers by publication support",
            os.path.join(figdir, "dotplot_celltype_markers_publication_support"),
            grouped=True,
            y_label="Cell type",
            row_column="cell_type",
            row_order_column="cell_type_order",
            size_column="n_pub",
            size_label="n.pub",
        )
    print(f"  [plots] {tag}: rendered reference-free figures -> {figdir}")


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default=root / "results", type=Path)
    parser.add_argument("--dataset", action="append", default=[])
    arguments = parser.parse_args()
    requested = set(arguments.dataset)
    result_dirs = [
        path
        for path in sorted(arguments.results_dir.iterdir())
        if path.is_dir()
        and (path / FIGDATA_DIR / "cells.csv").is_file()
        and (not requested or path.name in requested)
    ]
    if requested - {path.name for path in result_dirs}:
        missing = ", ".join(sorted(requested - {path.name for path in result_dirs}))
        raise FileNotFoundError(f"result datasets not found: {missing}")
    for result_dir in result_dirs:
        plot_dataset(result_dir.name, result_dir)
