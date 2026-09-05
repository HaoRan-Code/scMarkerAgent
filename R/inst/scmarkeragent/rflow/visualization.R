# =============================================================================
# visualization.R -- R-native rendering of the reference-free production figures.
#
# Renders the file set every arm renders, from the same figure_data inputs, with the
# same categories, category orders, palettes, dimensions and axis semantics. Bytes are
# not the bar -- two graphics stacks never agree on bytes -- so
# what this file guards is that a reader sees the SAME figure, not the same file.
#
# Reads ONLY figure_data/*.csv and config/plot_params.json. Never touches run caches,
# never re-decides an annotation.
#
#   1 umap_cluster      UMAP colored by cluster id, labels on data
#   2 umap_celltype     UMAP colored by primary annotation; separate legend panel,
#                       plus the stitched '<stem>_labeled.png'
#   3 dotplot_celltype_markers                  identity markers, one row per cluster
#   4 dotplot_celltype_markers_by_celltype      same markers, clusters pooled by label
#   5 dotplot_celltype_markers_publication_support   dot area = curated n_pub (log)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

.vis_script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      value <- frame$ofile
      if (is.null(value)) "" else as.character(value)[1]
    },
    ""
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files)) {
    return(dirname(normalizePath(frame_files[length(frame_files)])))
  }
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("^--file=", args, value = TRUE)
  if (length(match)) {
    return(dirname(normalizePath(sub("^--file=", "", match[1]))))
  }
  normalizePath(getwd())
}
.VIS_SELF <- .vis_script_dir()

if (!exists("scma_pp_num", mode = "function")) {
  source(file.path(.VIS_SELF, "plot_params.R"))
}

VIS_FIGDATA_DIR <- "figure_data"
VIS_FIGURES_DIR <- "figures"

# ---- eviden palette defaults ----------------------------------------------------
.EVIDEN_DEFAULT <- c(
  "#a4cde1", "#67a4cc", "#277fb8", "#549da3", "#96cb8f", "#8bc96d", "#4dae47",
  "#5c9e43", "#b79973", "#f38989", "#ec5051", "#e32427", "#ef6a45", "#f9b769",
  "#f9a341", "#f48521", "#ee8e46", "#d4a6a8", "#af93c4", "#8660a8", "#815e99",
  "#c6b598", "#f6f28f", "#d4a55b", "#b05a28"
)

# The tab20 / tab20b / tab20c qualitative sets, spelled out because the extension order
# is part of the palette contract every arm keeps.
.TAB20 <- c(
  "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728",
  "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2",
  "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
)
.TAB20B <- c(
  "#393b79", "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b",
  "#cedb9c", "#8c6d31", "#bd9e39", "#e7ba52", "#e7cb94", "#843c39", "#ad494a",
  "#d6616b", "#e7969c", "#7b4173", "#a55194", "#ce6dbd", "#de9ed6"
)
.TAB20C <- c(
  "#3182bd", "#6baed6", "#9ecae1", "#c6dbef", "#e6550d", "#fd8d3c", "#fdae6b",
  "#fdd0a2", "#31a354", "#74c476", "#a1d99b", "#c7e9c0", "#756bb1", "#9e9ac8",
  "#bcbddc", "#dadaeb", "#636363", "#969696", "#bdbdbd", "#d9d9d9"
)

# ---- shared figure params (EDIT config/plot_params.json, NOT this file) ----------
.VIS_FONT <- scma_pp_text("global", "font_family", "Arial")
.VIS_BASE_FS <- scma_pp_num("global", "base_font_size", 8, lo = 1)
.VIS_DPI <- scma_pp_num("global", "dpi", 300, lo = 1)
.VIS_FORMATS <- tolower(scma_pp_str_list("global", "formats", c("pdf", "png")))
.VIS_BG <- scma_pp_text("global", "figure_background", "white")

.VIS_EVIDEN_BASE <- scma_pp_str_list("palette", "eviden_base", .EVIDEN_DEFAULT)
.VIS_GRADIENT <- scma_pp_str_list(
  "dotplot", "gradient_colors", c("#fee0d2", "#fc9272", "#de2d26")
)

.PT_PER_MM <- 72.27 / 25.4

#' n distinct colors: eviden base first, then tab20/tab20b/tab20c hues, then HSV
#' fallback -- the extension order every arm uses.
scma_eviden_palette <- function(n) {
  if (n <= length(.VIS_EVIDEN_BASE)) {
    return(.VIS_EVIDEN_BASE[seq_len(n)])
  }
  palette <- .VIS_EVIDEN_BASE
  for (extension in list(.TAB20, .TAB20B, .TAB20C)) {
    palette <- c(palette, setdiff(extension, palette))
    if (length(palette) >= n) break
  }
  if (length(palette) < n) {
    hues <- seq(0, 1, length.out = n - length(palette) + 1L)[seq_len(n - length(palette))]
    palette <- c(palette, grDevices::hsv(hues, 1, 1))
  }
  palette[seq_len(n)]
}

#' Dynamic point size: area (pt^2) from the cell count, then clamped.
.vis_point_area <- function(n) {
  fixed <- scma_pp_num("umap", "point_size", 0, lo = 0)
  if (length(fixed) && fixed > 0) {
    return(as.numeric(fixed))
  }
  max(1.5, min(14.0, 90000.0 / max(n, 1)))
}

#' A scatter point is specified as an AREA in pt^2; ggplot's size is closer to a
#' DIAMETER in mm.
.vis_area_to_size_mm <- function(area_pt2) {
  2 * sqrt(pmax(as.numeric(area_pt2), 0) / pi) / .PT_PER_MM
}

#' Dot AREA (pt^2) for a percentage, following ggplot2's area_pal -- endpoints give
#' smin at 0 and smax at 100.
.vis_dot_area <- function(pct, smin, smax, cap) {
  v <- pmin(pmax(as.numeric(pct) / 100, 0), 1)
  pmin((sqrt(smin) + (sqrt(smax) - sqrt(smin)) * sqrt(v))^2, cap)
}

.vis_theme <- function(base_size = .VIS_BASE_FS) {
  theme_classic(base_size = base_size, base_family = .VIS_FONT) +
    theme(
      axis.line = element_line(linewidth = 0.6 / 2.14),
      axis.ticks = element_line(linewidth = 0.6 / 2.14),
      plot.title = element_text(face = "bold", size = base_size),
      axis.title = element_text(size = base_size),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.background = element_rect(fill = .VIS_BG, colour = NA),
      panel.background = element_rect(fill = .VIS_BG, colour = NA)
    )
}

#' Save one figure in every configured format at an exact device size in mm.
.vis_save <- function(draw, stem, width_mm, height_mm) {
  for (format in .VIS_FORMATS) {
    path <- sprintf("%s.%s", stem, format)
    if (format == "png") {
      grDevices::png(
        path,
        width = width_mm, height = height_mm, units = "mm",
        res = .VIS_DPI, type = "cairo", bg = .VIS_BG
      )
    } else if (format == "pdf") {
      grDevices::cairo_pdf(
        path,
        width = width_mm / 25.4, height = height_mm / 25.4, bg = .VIS_BG
      )
    } else {
      next
    }
    draw()
    grDevices::dev.off()
  }
}

# ------------------------------------------------------------------ UMAP panels ----

#' Shared SQUARE data limits for ALL UMAP panels of a dataset, computed once from every
#' cell -- panels showing subsets sit in the identical coordinate frame.
scma_umap_limits <- function(cells) {
  x <- as.numeric(cells$umap_x)
  y <- as.numeric(cells$umap_y)
  cx <- (min(x) + max(x)) / 2
  cy <- (min(y) + max(y)) / 2
  half <- 0.5 * max(max(x) - min(x), max(y) - min(y)) * 1.05
  if (!is.finite(half) || half == 0) half <- 1
  list(xlim = c(cx - half, cx + half), ylim = c(cy - half, cy + half))
}

#' Force the ggplot PANEL (the data area) to an exact square in mm -- the comparability
#' guarantee across panels and datasets.
.vis_fixed_panel <- function(plot, panel_mm) {
  gt <- ggplot2::ggplotGrob(plot)
  panel <- gt$layout[gt$layout$name == "panel", ]
  gt$widths[[panel$l]] <- grid::unit(panel_mm, "mm")
  gt$heights[[panel$t]] <- grid::unit(panel_mm, "mm")
  gt
}

.vis_umap_panel <- function(cells, categories, color_of, title, stem,
                            on_data_labels, order, xlim, ylim) {
  point_size <- .vis_area_to_size_mm(.vis_point_area(nrow(cells)))
  alpha <- scma_pp_num("umap", "point_alpha", 1.0, lo = 0, hi = 1)
  rasterize <- scma_pp_boolean("umap", "rasterize", TRUE)
  label_fs <- scma_pp_num("umap", "label_font_size", 7, lo = 1)
  panel_mm <- scma_pp_num("umap", "panel_mm", 55, lo = 10)
  m_l <- scma_pp_num("umap", "margin_left_mm", 8, lo = 0)
  m_r <- scma_pp_num("umap", "margin_right_mm", 3, lo = 0)
  m_b <- scma_pp_num("umap", "margin_bottom_mm", 7, lo = 0)
  m_t <- scma_pp_num("umap", "margin_top_mm", 9, lo = 0)
  width_mm <- panel_mm + m_l + m_r
  height_mm <- panel_mm + m_b + m_t

  frame <- data.table(
    umap_x = as.numeric(cells$umap_x),
    umap_y = as.numeric(cells$umap_y),
    category = factor(as.character(categories), levels = order)
  )
  points <- geom_point(
    aes(colour = category),
    size = point_size, alpha = alpha, stroke = 0, shape = 16
  )
  if (isTRUE(rasterize) && requireNamespace("ggrastr", quietly = TRUE)) {
    points <- ggrastr::rasterise(points, dpi = .VIS_DPI)
  }
  plot <- ggplot(frame, aes(umap_x, umap_y)) +
    points +
    scale_colour_manual(values = color_of, guide = "none") +
    coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on") +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    .vis_theme() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.margin = grid::unit(c(1, 1, 1, 1), "mm")
    )
  if (isTRUE(on_data_labels)) {
    medians <- frame[
      !is.na(category), .(
        umap_x = stats::median(umap_x), umap_y = stats::median(umap_y)
      ),
      by = category
    ]
    plot <- plot + geom_label(
      data = medians,
      aes(umap_x, umap_y, label = category),
      inherit.aes = FALSE,
      size = label_fs / .pt,
      fontface = "bold",
      colour = "black",
      fill = scales::alpha("white", 0.65),
      linewidth = 0,
      label.padding = grid::unit(0.1, "lines")
    )
  }
  # The gtable is built INSIDE the draw closure: ggplotGrob measures text on the
  # currently open device, and building it outside would measure against the null pdf()
  # device -- which lacks Arial and warns once per string.
  .vis_save(
    function() {
      grid::grid.newpage()
      grid::grid.draw(.vis_fixed_panel(plot, panel_mm))
    },
    stem, width_mm, height_mm
  )
}

#' The category legend as its OWN standalone panel, sized to its entries. Drawn from
#' data rather than extracted from a plot, so it cannot silently vanish when a ggplot
#' internal changes.
.vis_umap_legend <- function(order, color_of, legend_title, stem) {
  marker_size <- scma_pp_num("umap", "legend_marker_size", 4, lo = 0.1)
  ncol <- as.integer(scma_pp_num("umap", "legend_ncol", 1, lo = 1))
  n_rows <- ceiling(length(order) / max(ncol, 1))
  longest <- max(c(nchar(order), 8))
  width_mm <- 10 + ncol * (6 + longest * 1.5)
  height_mm <- max(20, 6 + n_rows * 4.4)

  entries <- data.table(
    label = factor(order, levels = order),
    column = (seq_along(order) - 1L) %/% n_rows,
    row = (seq_along(order) - 1L) %% n_rows
  )
  plot <- ggplot(entries) +
    geom_point(
      aes(x = column, y = -row, colour = label),
      # legend_marker_size is a marker DIAMETER in points.
      size = marker_size / .PT_PER_MM,
      shape = 16
    ) +
    geom_text(
      aes(x = column + 0.08, y = -row, label = label),
      hjust = 0, size = (.VIS_BASE_FS - 2) / .pt, family = .VIS_FONT
    ) +
    scale_colour_manual(values = color_of, guide = "none") +
    scale_x_continuous(limits = c(-0.15, max(entries$column) + 1)) +
    scale_y_continuous(limits = c(-(n_rows - 0.25), 1.25)) +
    labs(title = legend_title) +
    theme_void(base_family = .VIS_FONT) +
    theme(
      plot.title = element_text(size = .VIS_BASE_FS - 1, hjust = 0),
      plot.margin = grid::unit(c(2, 2, 2, 3), "mm"),
      plot.background = element_rect(fill = .VIS_BG, colour = NA)
    )
  .vis_save(function() print(plot), stem, width_mm, height_mm)
}

#' Horizontally stitch the fixed data panel + its legend into one PNG, panel left,
#' legend right, vertically centered -- the delivered `<stem>_labeled.png`.
.vis_stitch_umap <- function(panel_png, legend_png, out_png) {
  if (!file.exists(panel_png) || !file.exists(legend_png)) {
    return(invisible(NULL))
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    return(invisible(NULL))
  }
  pad_to <- function(image, height) {
    if (dim(image)[1] == height) {
      return(image)
    }
    top <- (height - dim(image)[1]) %/% 2
    canvas <- array(1, dim = c(height, dim(image)[2], dim(image)[3]))
    canvas[top + seq_len(dim(image)[1]), , ] <- image
    canvas
  }
  normalise <- function(image) {
    if (length(dim(image)) == 2L) {
      image <- array(rep(image, 3), dim = c(dim(image), 3L))
    }
    image[, , 1:3, drop = FALSE]
  }
  a <- normalise(png::readPNG(panel_png))
  b <- normalise(png::readPNG(legend_png))
  height <- max(dim(a)[1], dim(b)[1])
  a <- pad_to(a, height)
  b <- pad_to(b, height)
  stitched <- array(1, dim = c(height, dim(a)[2] + dim(b)[2], 3L))
  stitched[, seq_len(dim(a)[2]), ] <- a
  stitched[, dim(a)[2] + seq_len(dim(b)[2]), ] <- b
  png::writePNG(stitched, out_png, dpi = .VIS_DPI)
  invisible(out_png)
}

.vis_umap <- function(cells, categories, color_of, title, stem,
                      on_data_labels = FALSE, legend_title = NULL,
                      order = NULL, xlim = NULL, ylim = NULL) {
  if (is.null(order)) {
    order <- sort(unique(as.character(categories)), method = "radix")
  }
  .vis_umap_panel(
    cells, categories, color_of, title, stem,
    on_data_labels = on_data_labels, order = order, xlim = xlim, ylim = ylim
  )
  if (!on_data_labels && scma_pp_boolean("umap", "separate_legend", TRUE)) {
    .vis_umap_legend(order, color_of, legend_title, paste0(stem, "_legend"))
    .vis_stitch_umap(
      paste0(stem, ".png"), paste0(stem, "_legend.png"), paste0(stem, "_labeled.png")
    )
  }
}

# --------------------------------------------------------------- dotplot panel ----

#' Representative fill for the size-legend keys: midpoint of the fill gradient.
.vis_gradient_mid <- function() {
  override <- scma_pp_text("dotplot", "legend_dot_fill", "")
  if (nzchar(override)) {
    return(override)
  }
  ramp <- grDevices::colorRampPalette(.VIS_GRADIENT)(101)
  ramp[51]
}

#' Unified dotplot renderer -- the counterpart of visualization._dotplot. When `grouped`,
#' genes are laid out block-diagonal with separators and per-block labels; the size
#' channel is pct_exp, or a publication count placed on a log scale.
.vis_dotplot <- function(long, title, stem, grouped = FALSE,
                         y_label = "Cluster: cell type",
                         row_column = "cluster_celltype",
                         row_order_column = "cluster_order",
                         size_column = "pct_exp",
                         size_label = "pct.exp") {
  long <- as.data.table(long)
  if (is.null(long) || !nrow(long)) {
    cat(sprintf("  [plots] dotplot '%s': no data; skip\n", title))
    return(invisible(NULL))
  }
  size_range <- scma_pp_num_list("dotplot", "dot_size_range", c(1.0, 55.0), n = 2, lo = 0)
  smin <- size_range[1]
  smax <- size_range[2]
  cap <- scma_pp_num("dotplot", "dot_size_max_cap", 55.0, lo = 0)
  stroke <- scma_pp_num("dotplot", "dot_stroke", 0.9, lo = 0)
  tick_fs <- scma_pp_num("dotplot", "tick_font_size", 6, lo = 1)
  x_angle <- scma_pp_num("dotplot", "x_axis_angle", 90)
  grid_on <- scma_pp_boolean("dotplot", "panel_grid", TRUE)
  grid_alpha <- scma_pp_num("dotplot", "panel_grid_alpha", 0.3, lo = 0, hi = 1)
  max_genes <- as.integer(scma_pp_num("dotplot", "max_genes", 60, lo = 0))
  w_per <- scma_pp_num("dotplot", "width_mm_per_gene", 3.4)
  w_base <- scma_pp_num("dotplot", "width_mm_base", 55)
  w_min <- scma_pp_num("dotplot", "width_mm_min", 120)
  h_per <- scma_pp_num("dotplot", "height_mm_per_cluster", 6.0)
  h_base <- scma_pp_num("dotplot", "height_mm_base", 32)
  h_min <- scma_pp_num("dotplot", "height_mm_min", 60)
  block_gap <- scma_pp_num("dotplot", "block_gap", 0.9, lo = 0)
  sep_on <- scma_pp_boolean("dotplot", "block_separator", TRUE)
  sep_alpha <- scma_pp_num("dotplot", "block_separator_alpha", 0.55, lo = 0, hi = 1)
  lab_on <- scma_pp_boolean("dotplot", "show_block_labels", TRUE)
  lab_fs <- scma_pp_num("dotplot", "block_label_font_size", 5.5, lo = 1)
  lab_ang <- scma_pp_num("dotplot", "block_label_angle", 30)
  lab_max <- as.integer(scma_pp_num("dotplot", "block_label_max_chars", 22, lo = 1))

  has_group <- isTRUE(grouped) && ("gene_group" %in% names(long)) &&
    any(nzchar(as.character(long$gene_group)))

  slot_column <- if ("marker_slot" %in% names(long)) "marker_slot" else "gene"
  slot_rows <- unique(
    long[base::order(long$gene_order)],
    by = slot_column
  )
  slots <- as.character(slot_rows[[slot_column]])
  slot_labels <- setNames(as.character(slot_rows$gene), slots)
  if (max_genes > 0 && length(slots) > max_genes) {
    slots <- slots[seq_len(max_genes)]
    long <- long[long[[slot_column]] %in% slots]
  }
  row_frame <- unique(
    long[base::order(long[[row_order_column]])],
    by = row_column
  )
  rows <- as.character(row_frame[[row_column]])
  row_index <- setNames(seq_along(rows) - 1L, rows)

  size_values <- as.numeric(long[[size_column]])
  size_max <- suppressWarnings(max(size_values, na.rm = TRUE))
  if (!is.finite(size_max) || size_max <= 0) size_max <- 1
  if (size_column == "pct_exp") {
    size_fraction <- function(value) as.numeric(value) / size_max
    size_keys <- c(0, 25, 50, 75, 100)
  } else {
    size_fraction <- function(value) {
      counts <- pmax(as.numeric(value), 0)
      log1p(counts) / log1p(size_max)
    }
    size_keys <- sort(unique(vapply(
      c(0, 0.25, 0.5, 0.75, 1),
      function(f) as.integer(round(expm1(f * log1p(size_max)))),
      0L
    )))
  }

  block_of <- setNames(rep(0L, length(slots)), slots)
  label_of_slot <- setNames(rep("", length(slots)), slots)
  if (has_group) {
    groups <- unique(long, by = slot_column)
    block_of[as.character(groups[[slot_column]])] <-
      as.integer(groups$gene_group_order)
    label_of_slot[as.character(groups[[slot_column]])] <-
      as.character(groups$gene_group)
  }
  x_of <- setNames(
    (seq_along(slots) - 1L) + if (has_group) block_gap * block_of[slots] else 0,
    slots
  )
  x_max <- if (length(x_of)) max(x_of) else 0
  fill_limits <- range(as.numeric(long$avg_exp_scaled))
  if (fill_limits[1] == fill_limits[2]) {
    fill_limits <- fill_limits + c(-1e-6, 1e-6)
  }

  height_mm <- max(h_min, length(rows) * h_per + h_base)
  width_mm <- max(w_min, (x_max + 1) * w_per + w_base)

  frame <- data.table(
    x = unname(x_of[as.character(long[[slot_column]])]),
    y = unname(row_index[as.character(long[[row_column]])]),
    fill_value = as.numeric(long$avg_exp_scaled),
    area = .vis_dot_area(100 * size_fraction(size_values), smin, smax, cap)
  )
  frame[, size_mm := .vis_area_to_size_mm(area)]

  # The legend keys are FIXED percentages / counts, not data quantiles. An identity
  # scale drops breaks outside its limits, and the data rarely spans the whole range,
  # so the limits are widened to hold every key.
  key_sizes <- .vis_area_to_size_mm(
    .vis_dot_area(100 * size_fraction(size_keys), smin, smax, cap)
  )
  size_limits <- range(c(frame$size_mm, key_sizes))

  plot <- ggplot(frame, aes(x, y)) +
    geom_point(
      aes(fill = fill_value, size = size_mm),
      shape = 21, colour = "black", stroke = stroke / 2
    ) +
    scale_fill_gradientn(
      colours = .VIS_GRADIENT, limits = fill_limits, name = "avg.exp.scaled",
      # The bar keeps a fixed physical size so the whole legend column fits beside
      # even the shortest (h_min) panel instead of being clipped.
      guide = guide_colorbar(
        order = 1,
        barheight = grid::unit(
          scma_pp_num("dotplot", "colorbar_length_mm", 18, lo = 1), "mm"
        ),
        barwidth = grid::unit(
          scma_pp_num("dotplot", "colorbar_thickness_mm", 3.0, lo = 0.1), "mm"
        )
      )
    ) +
    scale_size_identity(
      guide = guide_legend(
        order = 2,
        keyheight = grid::unit(4.2, "mm"),
        keywidth = grid::unit(4.2, "mm"),
        override.aes = list(
          fill = .vis_gradient_mid(), colour = "black", stroke = stroke / 2
        )
      ),
      breaks = key_sizes,
      labels = as.character(size_keys),
      limits = size_limits,
      name = size_label
    ) +
    scale_x_continuous(
      breaks = unname(x_of[slots]),
      labels = unname(slot_labels[slots]),
      limits = c(-0.6, x_max + 0.6),
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      breaks = seq_along(rows) - 1L,
      labels = rows,
      limits = c(length(rows) - 0.4, -0.6),
      expand = c(0, 0)
    ) +
    labs(x = "Markers", y = y_label, title = title) +
    .vis_theme()

  # When per-block cell-type labels are drawn above the panel, push the title up by
  # their worst-case vertical rise so it never overlaps them -- the same reservation
  # every arm makes through its title pad (in points).
  title_pad_pt <- 6
  if (has_group && isTRUE(lab_on)) {
    title_pad_pt <- title_pad_pt +
      abs(sin(lab_ang * pi / 180)) * lab_max * lab_fs * 0.6 + lab_fs
  }
  plot <- plot +
    theme(
      axis.text.x = element_text(
        angle = x_angle, hjust = if (x_angle == 90) 1 else 0.5,
        vjust = if (x_angle == 90) 0.5 else 1, size = tick_fs
      ),
      axis.text.y = element_text(size = tick_fs),
      legend.title = element_text(size = .VIS_BASE_FS - 1),
      legend.text = element_text(size = tick_fs),
      legend.position = "right",
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = grid::unit(2, "mm"),
      legend.box.margin = margin(0, 0, 0, 1),
      plot.title = element_text(
        face = "bold", size = .VIS_BASE_FS,
        margin = margin(b = title_pad_pt, unit = "pt")
      ),
      panel.grid.major = if (grid_on) {
        element_line(colour = scales::alpha("grey", grid_alpha), linewidth = 0.3 / 2.14)
      } else {
        element_blank()
      }
    )

  if (has_group) {
    blocks <- split(unname(x_of[slots]), block_of[slots])
    ordered <- sort(as.integer(names(blocks)))
    if (isTRUE(sep_on) && length(ordered) > 1) {
      separators <- vapply(
        seq_len(length(ordered) - 1L),
        function(i) {
          (max(blocks[[as.character(ordered[i])]]) +
            min(blocks[[as.character(ordered[i + 1L])]])) / 2
        },
        0
      )
      plot <- plot + geom_vline(
        xintercept = separators,
        colour = scales::alpha("grey", sep_alpha), linewidth = 0.4 / 2.14
      )
    }
    if (isTRUE(lab_on)) {
      label_frame <- rbindlist(lapply(ordered, function(block) {
        xs <- blocks[[as.character(block)]]
        label <- label_of_slot[slots][block_of[slots] == block][1]
        if (nchar(label) > lab_max) {
          label <- paste0(substr(label, 1, lab_max - 1L), "\u2026")
        }
        data.table(x = (min(xs) + max(xs)) / 2, label = label)
      }))
      plot <- plot +
        geom_text(
          data = label_frame,
          aes(x = x, y = -0.3, label = label),
          inherit.aes = FALSE,
          angle = lab_ang,
          size = lab_fs / .pt,
          hjust = if (lab_ang != 0) 0 else 0.5,
          vjust = 0,
          family = .VIS_FONT
        ) +
        coord_cartesian(clip = "off")
    }
  }

  .vis_save(function() print(plot), stem, width_mm, height_mm)
}

# ------------------------------------------------------------------- driver ----

#' Render reference-free figures for one completed production run, reading only the
#' figure_data bundle.
scma_plot_dataset <- function(tag, outdir, fig_subdir = VIS_FIGURES_DIR) {
  bundle <- file.path(outdir, VIS_FIGDATA_DIR)
  figdir <- if (nzchar(fig_subdir)) file.path(outdir, fig_subdir) else outdir
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

  read_bundle <- function(name, required = TRUE) {
    path <- file.path(bundle, name)
    if (!file.exists(path)) {
      if (required) stop(sprintf("missing figure data: %s", path), call. = FALSE)
      return(NULL)
    }
    fread(path, colClasses = list(character = intersect(
      c("cluster", "cluster_celltype", "cell_type", "gene_group", "marker_slot", "gene"),
      names(fread(path, nrows = 0))
    )), na.strings = NULL, encoding = "UTF-8")
  }

  cells <- read_bundle("cells.csv")
  clustermap <- read_bundle("clustermap.csv")
  identity_markers <- read_bundle("dotplot_celltype_markers.csv")
  celltype_markers <- read_bundle(
    "dotplot_celltype_markers_by_celltype.csv",
    required = FALSE
  )

  umap_ok <- ("umap_x" %in% names(cells)) &&
    any(is.finite(suppressWarnings(as.numeric(cells$umap_x))))
  if (!umap_ok) {
    cat(sprintf("  [plots] %s: no UMAP coords in bundle; skipping UMAP panels\n", tag))
  }

  order_map <- as.data.table(clustermap)[order(cluster_order)]
  clusters <- as.character(order_map$cluster)
  cc_of <- setNames(as.character(order_map$cluster_celltype), clusters)
  palette <- scma_eviden_palette(length(clusters))
  color_by_cluster <- setNames(palette, clusters)
  color_by_cc <- setNames(palette, unname(cc_of[clusters]))
  cc_order <- unname(cc_of[clusters])

  if (umap_ok) {
    limits <- scma_umap_limits(cells)
    .vis_umap(
      cells, as.character(cells$cluster), color_by_cluster, "Clusters",
      file.path(figdir, "umap_cluster"),
      on_data_labels = TRUE, order = clusters,
      xlim = limits$xlim, ylim = limits$ylim
    )
    cc_values <- unname(cc_of[as.character(cells$cluster)])
    .vis_umap(
      cells, cc_values, color_by_cc, "Primary annotation",
      file.path(figdir, "umap_celltype"),
      on_data_labels = FALSE, legend_title = "Cluster: primary annotation",
      order = cc_order, xlim = limits$xlim, ylim = limits$ylim
    )
  }
  .vis_dotplot(
    identity_markers,
    "LLM-validated cell-type identity markers",
    file.path(figdir, "dotplot_celltype_markers"),
    grouped = TRUE,
    y_label = "Cluster: primary annotation"
  )
  if (!is.null(celltype_markers)) {
    .vis_dotplot(
      celltype_markers,
      "Cell-type identity markers",
      file.path(figdir, "dotplot_celltype_markers_by_celltype"),
      grouped = TRUE,
      y_label = "Cell type",
      row_column = "cell_type",
      row_order_column = "cell_type_order"
    )
    .vis_dotplot(
      celltype_markers,
      "Cell-type identity markers by publication support",
      file.path(figdir, "dotplot_celltype_markers_publication_support"),
      grouped = TRUE,
      y_label = "Cell type",
      row_column = "cell_type",
      row_order_column = "cell_type_order",
      size_column = "n_pub",
      size_label = "n.pub"
    )
  }
  cat(sprintf("  [plots] %s: rendered reference-free figures -> %s\n", tag, figdir))
  invisible(figdir)
}
