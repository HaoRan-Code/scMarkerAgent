start_time <- Sys.time()
suppressPackageStartupMessages({
  library(UCell)
  library(BiocParallel)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(scales)  # for alpha() transparency in plots
  library(showtext)
  library(ggrastr)
})

# ---- UTF-8 / font configuration for ggplot2 ----
# Ensure that UTF-8 text (e.g. Greek letters like α, β) is rendered correctly in plots
font_add("Arial", "/usr/share/fonts/msttcore/arial.ttf")
showtext_auto()
showtext_opts(dpi = 300)

# ---- argument parsing (no extra deps) ----
.defaults <- list(
  species = "Human",
  tissue = "brain",
  disease_name = "glioblastoma",
  cell_name = "malignant cell",
  polarity = "negative",
  n_threads = "4",
  outdir = "./results/cell_score_example",
  input_seurat = "./data/example_glioblastoma_input.rds",
  marker_rds = "./data/mydata.RDS",
  type_in_metadata = "",
  target_cells = "",
  n_variable_features = "2000",
  dims = "30",
  k_param = "20",
  cluster_resolution = "0.5",
  random_seed = "1234"
)
parse_cli_args <- function(args, defaults) {
  opts <- defaults
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (grepl("^--", a)) {
      if (grepl("=", a, fixed = TRUE)) {
        kv <- sub("^--", "", a)
        key <- sub("=.*$", "", kv)
        val <- sub("^[^=]*=", "", kv)
        if (key %in% names(opts)) opts[[key]] <- val
        i <- i + 1
      } else {
        key <- sub("^--", "", a)
        if (i < length(args) && !grepl("^--", args[i + 1])) {
          val <- args[i + 1]
          if (key %in% names(opts)) opts[[key]] <- val
          i <- i + 2
        } else {
          i <- i + 1
        }
      }
    } else {
      i <- i + 1
    }
  }
  opts
}

args <- commandArgs(trailingOnly = TRUE)
opt <- parse_cli_args(args, .defaults)
opt$n_threads <- suppressWarnings(as.integer(opt$n_threads))
if (is.na(opt$n_threads) || opt$n_threads < 1) opt$n_threads <- 4

opt$n_variable_features <- as.integer(opt$n_variable_features)
opt$dims <- as.integer(opt$dims)
opt$k_param <- as.integer(opt$k_param)
opt$cluster_resolution <- as.numeric(opt$cluster_resolution)
opt$random_seed <- as.integer(opt$random_seed)

if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---- load inputs ----
example_ucell_input <- readRDS(opt$input_seurat)

# ---- preprocessing ----
message("Running standard Seurat preprocessing...")

assay <- "RNA"
DefaultAssay(example_ucell_input) <- assay

message("Normalizing data...")
example_ucell_input <- NormalizeData(example_ucell_input, assay = assay, verbose = FALSE)

message("Finding variable features...")
example_ucell_input <- FindVariableFeatures(
  example_ucell_input, 
  assay = assay, 
  selection.method = "vst", 
  nfeatures = opt$n_variable_features, 
  verbose = FALSE
)

message("Scaling data...")
example_ucell_input <- ScaleData(
  example_ucell_input, 
  assay = assay, 
  features = VariableFeatures(example_ucell_input), 
  verbose = FALSE
)

message("Running PCA...")
set.seed(opt$random_seed)
example_ucell_input <- RunPCA(
  example_ucell_input, 
  assay = assay, 
  features = VariableFeatures(example_ucell_input), 
  seed.use = opt$random_seed, 
  verbose = FALSE
)

message("Building neighborhood graph...")
example_ucell_input <- FindNeighbors(
  example_ucell_input, 
  dims = 1:opt$dims, 
  k.param = opt$k_param, 
  verbose = FALSE
)

message("Finding clusters...")
example_ucell_input <- FindClusters(
  example_ucell_input, 
  resolution = opt$cluster_resolution, 
  random.seed = opt$random_seed, 
  verbose = FALSE
)

message("Running UMAP...")
example_ucell_input <- RunUMAP(
  example_ucell_input, 
  dims = 1:opt$dims, 
  seed.use = opt$random_seed, 
  verbose = FALSE
)

mydata <- readRDS(opt$marker_rds)

required_cols <- c("species", "top_level", "disease_type_do", "cell_name_standardized", "marker_polarity", "gene_symbol")
missing_cols <- setdiff(required_cols, colnames(mydata))
if (length(missing_cols) > 0) {
  stop(sprintf("marker_rds is missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

# ---- build signature ----
cell_signatures_vec <- mydata %>%
  filter(
    species == opt$species,
    top_level == opt$tissue,
    disease_type_do == opt$disease_name,
    cell_name_standardized == opt$cell_name,
    marker_polarity == opt$polarity
  ) %>%
  .$gene_symbol %>%
  unique()

if (length(cell_signatures_vec) == 0) {
  stop("No marker genes matched the provided filters; please check species/tissue/disease_name/cell_name/polarity.")
}

cell_signatures <- list(cell_signatures = cell_signatures_vec)

# ---- score with UCell ----
bp <- MulticoreParam(workers = opt$n_threads)
example_ucell_input <- AddModuleScore_UCell(
  example_ucell_input,
  chunk.size = 5000,
  features = cell_signatures,
  name = NULL,
  BPPARAM = bp,
  force.gc = TRUE
)
gc()

# ---- rename score column for clarity ----
if ("cell_signatures" %in% colnames(example_ucell_input@meta.data)) {
  colnames(example_ucell_input@meta.data)[
    colnames(example_ucell_input@meta.data) == "cell_signatures"
  ] <- "cell_signatures_score"
}

# ---- append UMAP coordinates to metadata (if available) ----
reduction_names <- tryCatch(names(example_ucell_input@reductions), error = function(e) character(0))
if ("umap" %in% reduction_names) {
  umap_coords <- Embeddings(example_ucell_input, reduction = "umap")
  umap_df <- as.data.frame(umap_coords[, 1:2, drop = FALSE])
  if (ncol(umap_df) >= 2) {
    colnames(umap_df)[1:2] <- c("umap_1", "umap_2")
  } else if (ncol(umap_df) == 1) {
    colnames(umap_df)[1] <- "umap_1"
  }
  umap_df <- umap_df[rownames(example_ucell_input@meta.data), , drop = FALSE]
  example_ucell_input@meta.data <- cbind(example_ucell_input@meta.data, umap_df)
}

# ---- save tables and object ----
meta_out <- file.path(opt$outdir, "meta_with_ucell_scores.csv")
write.csv(example_ucell_input@meta.data, meta_out, row.names = TRUE)

# ---- plotting: fixed and optional custom ----
metadata <- example_ucell_input@meta.data
has_umap <- all(c("umap_1", "umap_2") %in% colnames(metadata))

if (!has_umap) {
  message("UMAP coordinates not found; skipping plots.")
} else {
  # Fixed plot: color by cell_signatures (continuous)
  p_fixed <- ggplot(metadata, aes(x = umap_1, y = umap_2, color = cell_signatures_score)) +
    ggrastr::rasterise(geom_point(size = 0.3), dpi = 300) +
    scale_color_viridis_c() +
    coord_equal() +
    theme_minimal() +
    theme(text = element_text(family = "Arial")) +
    labs(color = "Scores")
  ggplot2::ggsave(
    filename = file.path(opt$outdir, "FeaturePlot_cell_signatures_umap.png"),
    plot = p_fixed, width = 7, height = 5, dpi = 300
  )
  showtext_auto(FALSE)
  ggplot2::ggsave(
    filename = file.path(opt$outdir, "FeaturePlot_cell_signatures_umap.pdf"),
    plot = p_fixed, width = 7, height = 5, device = cairo_pdf
  )
  showtext_auto(TRUE)

  # Compute filtered column names (<100 unique values)
  colnames_filtered <- metadata %>%
    dplyr::select(where(~ length(unique(.)) < 100)) %>%
    colnames()

  ti <- opt$type_in_metadata
  if (!is.null(ti) && nchar(ti) > 0 && ti %in% colnames_filtered) {
    values <- metadata[[ti]]
    # Ensure categorical for manual palettes
    if (!is.factor(values)) values <- factor(values)
    n_levels <- nlevels(values)
    p_custom <- ggplot(metadata, aes(x = umap_1, y = umap_2, color = values)) +
      ggrastr::rasterise(geom_point(size = 0.3), dpi = 300) +
      coord_equal() +
      theme_minimal() +
      theme(text = element_text(family = "Arial")) +
      labs(color = ti) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    if (n_levels == 2) {
      p_custom <- p_custom +
        scale_color_manual(values = c("#ED5051", "#287FB6"))
    } else if (n_levels >= 3) {
      palette_multi <- c(
        "#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d",
        "#4dae47","#5c9e43","#b79973","#f38989","#ec5051","#e32427",
        "#ef6a45","#f9b769","#f9a341","#f48521","#ee8e46","#d4a6a8",
        "#af93c4","#8660a8","#815e99","#c6b598","#f6f28f","#d4a55b",
        "#b05a28"
      )
      p_custom <- p_custom +
        scale_color_manual(values = rep(palette_multi, length.out = n_levels))
    } else {
      # Fallback: single level, use a constant color
      p_custom <- p_custom + scale_color_manual(values = "#277fb8")
    }
    # Safe filename for type_in_metadata
    ti_safe <- gsub("[^A-Za-z0-9_\\-]+", "_", ti)
    ggplot2::ggsave(
      filename = file.path(opt$outdir, "CustomPlot_umap.png"),
      plot = p_custom, width = 7, height = 5, dpi = 300
    )
    showtext_auto(FALSE)
    ggplot2::ggsave(
      filename = file.path(opt$outdir, "CustomPlot_umap.pdf"),
      plot = p_custom, width = 7, height = 5, device = cairo_pdf
    )
    showtext_auto(TRUE)

    # Optional: binary target vs other based on target_cells value
    target_val <- opt$target_cells
    if (!is.null(target_val) && nchar(target_val) > 0) {
      is_target <- as.character(metadata[[ti]]) == target_val
      # Report table like the example
      tab <- table(is_target)
      message(sprintf("Binary split for '%s' == '%s': FALSE=%s TRUE=%s", ti, target_val, tab["FALSE"], tab["TRUE"]))
      target_label <- target_val
      other_label <- "other cells"
      target_binary <- ifelse(is_target, target_label, other_label)
      metadata$target_binary <- factor(target_binary, levels = c(target_label, other_label))

      # UMAP scatter colored by target vs other
      p_target <- ggplot(metadata, aes(x = umap_1, y = umap_2, color = target_binary)) +
        ggrastr::rasterise(geom_point(size = 0.3), dpi = 300) +
        coord_equal() +
        theme_minimal() +
        theme(text = element_text(family = "Arial")) +
        scale_color_manual(values = setNames(c("#B71B1B", "#00468B"), c(target_label, other_label))) +
        labs(color = "Group") +
        guides(color = guide_legend(override.aes = list(size = 3)))
      tv_safe <- gsub("[^A-Za-z0-9_\\-]+", "_", target_val)
      ggplot2::ggsave(
        filename = file.path(opt$outdir, "TargetCellsPlot_umap.png"),
        plot = p_target, width = 7, height = 5, dpi = 300
      )
      showtext_auto(FALSE)
      ggplot2::ggsave(
        filename = file.path(opt$outdir, "TargetCellsPlot_umap.pdf"),
        plot = p_target, width = 7, height = 5, device = cairo_pdf
      )
      showtext_auto(TRUE)

      # Boxplot: cell_signatures by target vs other
      # Define theme colors (opaque) following demo_boxplot_theme_fixed.R style
      theme_colors <- setNames(
        c("#B71B1B", "#00468B"),
        c(target_label, other_label)
      )
      # Apply 80% opacity for fill only, keep lines opaque
      theme_colors_alpha <- scales::alpha(theme_colors, 0.8)

      p_box <- ggplot(
        metadata,
        aes(
          x = target_binary,
          y = cell_signatures_score,
          fill = target_binary,
          color = target_binary
        )
      ) +
        geom_boxplot(width = 0.5, outlier.size = 2, size = 1.5) +
        scale_fill_manual(values = theme_colors_alpha) +
        scale_color_manual(values = theme_colors) +
        theme_minimal() +
        theme(
          text = element_text(family = "Arial", size = 12),
          axis.title.x = element_blank(),
          legend.position = "none"
        ) +
        labs(
          y = "Scores"
        )
      # Statistical test (Wilcoxon rank-sum) and star annotation
      group_a <- metadata$cell_signatures_score[metadata$target_binary == target_label]
      group_b <- metadata$cell_signatures_score[metadata$target_binary == other_label]
      valid_a <- sum(!is.na(group_a)) > 0
      valid_b <- sum(!is.na(group_b)) > 0
      if (valid_a && valid_b) {
        wt <- tryCatch(
          wilcox.test(group_a, group_b, alternative = "two.sided"),
          error = function(e) NULL
        )
        if (!is.null(wt)) {
          pval <- wt$p.value
          stars <- if (pval < 1e-4) {
            "****"
          } else if (pval < 1e-3) {
            "***"
          } else if (pval < 1e-2) {
            "**"
          } else if (pval < 5e-2) {
            "*"
          } else {
            "ns"
          }
          message(sprintf("Wilcoxon rank-sum p-value (target vs other): %.3g (%s)", pval, stars))
          y_max <- suppressWarnings(max(metadata$cell_signatures_score, na.rm = TRUE))
          y_span <- if (is.finite(y_max)) abs(y_max) else 1
          y_pos <- y_max + 0.05 * (y_span + 1e-6)
          y_top <- y_pos + 0.03 * (y_span + 1e-6)
          p_box <- p_box +
            expand_limits(y = y_top) +
            annotate("segment", x = 1, xend = 2, y = y_pos, yend = y_pos) +
            annotate("segment", x = 1, xend = 1, y = y_pos, yend = y_pos - 0.02 * (y_span + 1e-6)) +
            annotate("segment", x = 2, xend = 2, y = y_pos, yend = y_pos - 0.02 * (y_span + 1e-6)) +
            annotate("text", x = 1.5, y = y_pos + 0.02 * (y_span + 1e-6), label = stars)
        }
      }
      ggplot2::ggsave(
        filename = file.path(opt$outdir, "Boxplot_cell_signatures_by_target.png"),
        plot = p_box, width = 6, height = 4, dpi = 300
      )
      showtext_auto(FALSE)
      ggplot2::ggsave(
        filename = file.path(opt$outdir, "Boxplot_cell_signatures_by_target.pdf"),
        plot = p_box, width = 6, height = 4, device = cairo_pdf
      )
      showtext_auto(TRUE)
    }
  } else if (!is.null(ti) && nchar(ti) > 0) {
    message(sprintf("type_in_metadata '%s' is not in colnames_filtered; skipping custom plot and only saving fixed plot.", ti))
  }
}

message(sprintf("UCell scoring completed. Outputs saved to: %s", normalizePath(opt$outdir)))
end_time <- Sys.time()
execution_time <- end_time - start_time
cat("Total execution time: ", execution_time, "\n")
