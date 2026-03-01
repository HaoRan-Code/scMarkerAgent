# ============================================================================
# Visualization Functions
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(ggrastr)
})

# Base color palette
.scmarkeragent_base_colors <- c(
  "#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d",
  "#4dae47","#5c9e43","#b79973","#f38989","#ec5051","#e32427",
  "#ef6a45","#f9b769","#f9a341","#f48521","#ee8e46","#d4a6a8",
  "#af93c4","#8660a8","#815e99","#c6b598","#f6f28f","#d4a55b",
  "#b05a28"
)

# Build a qualitative palette of length n using:
# 1) user-provided colors
# 2) RColorBrewer qualitative palettes to extend
# 3) fallback to distinct hues if still insufficient
.scmarkeragent_build_palette <- function(n) {
  if (n <= length(.scmarkeragent_base_colors)) {
    return(.scmarkeragent_base_colors[seq_len(n)])
  }
  pal <- .scmarkeragent_base_colors
  brewer_sets <- c("Paired","Set3","Dark2","Set2","Accent","Set1","Pastel1","Pastel2")
  for (set_name in brewer_sets) {
    if (set_name %in% rownames(RColorBrewer::brewer.pal.info)) {
      max_n <- RColorBrewer::brewer.pal.info[set_name, "maxcolors"]
      pal <- c(pal, RColorBrewer::brewer.pal(max_n, set_name))
      if (length(pal) >= n) break
    }
  }
  if (length(pal) < n) {
    # fallback: add distinct hues
    needed <- n - length(pal)
    pal <- c(pal, scales::hue_pal()(needed))
  }
  pal[seq_len(n)]
}

#' Generate cluster and annotation DimPlots
#' 
#' @param seurat_obj Seurat object
#' @param cluster_col Column name for cluster assignments
#' @param annotation_col Column name for cell type annotations
#' @param output_prefix Output file path prefix
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @return NULL (saves files to disk)
plot_cluster_and_annotation <- function(
  seurat_obj,
  cluster_col,
  annotation_col,
  output_prefix,
  width = 14,
  height = 6
) {
  
  # Sort clusters numerically if possible
  cluster_values <- as.character(seurat_obj@meta.data[[cluster_col]])
  cluster_numeric <- suppressWarnings(as.numeric(cluster_values))
  
  if (!any(is.na(cluster_numeric))) {
    cluster_levels <- as.character(sort(unique(cluster_numeric)))
    seurat_obj@meta.data[[cluster_col]] <- factor(
      cluster_values, 
      levels = cluster_levels
    )
  }
  
  # Determine palette sizes
  n_clusters <- length(unique(seurat_obj@meta.data[[cluster_col]]))
  n_annots <- length(unique(seurat_obj@meta.data[[annotation_col]]))
  cluster_cols <- .scmarkeragent_build_palette(n_clusters)
  annot_cols <- .scmarkeragent_build_palette(n_annots)
  
  # Cluster plot (with text labels)
  p1 <- DimPlot(
    seurat_obj, 
    group.by = cluster_col,
    cols = cluster_cols,
    label = TRUE, 
    repel = TRUE,
    raster = TRUE
  ) +
    ggtitle("Clusters") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.text = element_text(size = 10)
    )
  
  # Annotation plot (with text labels)
  p2 <- DimPlot(
    seurat_obj, 
    group.by = annotation_col,
    cols = annot_cols,
    label = TRUE, 
    repel = TRUE,
    raster = TRUE
  ) +
    ggtitle("Cell Type Annotation") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.text = element_text(size = 10)
    )
  
  # Combined plot (labeled version)
  p_combined <- p1 | p2
  
  # No-label version (legend only)
  p1_no_text <- DimPlot(
    seurat_obj, 
    group.by = cluster_col,
    cols = cluster_cols,
    label = FALSE, 
    repel = FALSE,
    raster = TRUE
  ) +
    ggtitle("Clusters") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.text = element_text(size = 10)
    )
  
  p2_no_text <- DimPlot(
    seurat_obj, 
    group.by = annotation_col,
    cols = annot_cols,
    label = FALSE, 
    repel = FALSE,
    raster = TRUE
  ) +
    ggtitle("Cell Type Annotation") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.text = element_text(size = 10)
    )
  
  p_combined_no_text <- p1_no_text | p2_no_text
  
  # Save PDF (labeled)
  showtext_auto(FALSE)
  ggsave(
    paste0(output_prefix, "_cluster_annotation.pdf"),
    p_combined, 
    width = width, 
    height = height,
    device = cairo_pdf
  )
  showtext_auto(TRUE)
  
  # Save PNG (labeled)
  ggsave(
    paste0(output_prefix, "_cluster_annotation.png"),
    p_combined, 
    width = width, 
    height = height,
    dpi = 300,
    bg = "white"
  )
  
  # Save PDF (no labels)
  showtext_auto(FALSE)
  ggsave(
    paste0(output_prefix, "_cluster_annotation_no_text.pdf"),
    p_combined_no_text, 
    width = width, 
    height = height,
    device = cairo_pdf
  )
  showtext_auto(TRUE)
  
  # Save PNG (no labels)
  ggsave(
    paste0(output_prefix, "_cluster_annotation_no_text.png"),
    p_combined_no_text, 
    width = width, 
    height = height,
    dpi = 300,
    bg = "white"
  )
  
  cat(sprintf("  ✓ Saved: %s_cluster_annotation.pdf/png (with text)\n", basename(output_prefix)))
  cat(sprintf("  ✓ Saved: %s_cluster_annotation_no_text.pdf/png (no text)\n", basename(output_prefix)))
}


#' Generate marker dot plot (LLM-validated key markers)
#' 
#' @param seurat_obj Seurat object
#' @param llm_results LLM decision results (must contain selected_celltype and key_markers_validated)
#' @param cluster_col Column name for cluster assignments
#' @param output_prefix Output file path prefix
#' @param max_markers_per_celltype Maximum number of markers to display per cell type
#' @return NULL (saves files to disk)
plot_marker_dotplot <- function(
  seurat_obj,
  llm_results,
  cluster_col,
  output_prefix,
  max_markers_per_celltype = 10
) {
  
  # Map cluster -> celltype -> markers
  celltype_markers <- list()
  cluster_to_celltype <- setNames(
    llm_results$selected_celltype,
    llm_results$cluster
  )
  
  # Assign each cell its annotated cell type
  seurat_obj@meta.data$temp_celltype_annotation <- 
    cluster_to_celltype[as.character(seurat_obj@meta.data[[cluster_col]])]
  
  # Collect key markers per cell type
  for (i in seq_len(nrow(llm_results))) {
    celltype <- llm_results$selected_celltype[i]
    markers_str <- llm_results$key_markers_validated[i]
    
    # Skip Unknown
    if (celltype == "Unknown" || is.na(celltype)) {
      next
    }
    
    # Parse markers (comma- or semicolon-separated)
    if (is.na(markers_str) || markers_str == "" || markers_str == "character(0)") {
      next
    }
    
    # Try multiple delimiters
    markers <- unlist(strsplit(markers_str, "[,;]"))
    markers <- trimws(markers)
    markers <- markers[markers != ""]
    
    if (length(markers) == 0) {
      next
    }
    
    # Cap markers per cell type
    if (length(markers) > max_markers_per_celltype) {
      markers <- head(markers, max_markers_per_celltype)
    }
    
    # Merge markers for the same cell type
    if (celltype %in% names(celltype_markers)) {
      celltype_markers[[celltype]] <- unique(c(celltype_markers[[celltype]], markers))
    } else {
      celltype_markers[[celltype]] <- markers
    }
  }
  
  # Check for valid markers
  if (length(celltype_markers) == 0) {
    cat("  ⚠ No valid markers found for DotPlot\n")
    return(invisible(NULL))
  }
  
  # Determine cell type order (based on cluster order in llm_results)
  ordered_celltypes <- unique(llm_results$selected_celltype)
  ordered_celltypes <- ordered_celltypes[!is.na(ordered_celltypes) & 
                                         ordered_celltypes != "Unknown" & 
                                         ordered_celltypes != ""]
  
  # Determine marker order (grouped by cell type order)
  ordered_markers <- c()
  for (ct in ordered_celltypes) {
    if (ct %in% names(celltype_markers)) {
      ct_markers <- celltype_markers[[ct]]
      ordered_markers <- c(ordered_markers, ct_markers)
    }
  }
  
  # Deduplicate while preserving first-occurrence order
  ordered_markers <- unique(ordered_markers)
  all_markers <- ordered_markers
  
  # Intersect with features present in the Seurat object (preserves order)
  available_features <- rownames(seurat_obj)
  valid_markers <- intersect(ordered_markers, available_features)
  
  if (length(valid_markers) == 0) {
    cat("  ⚠ No markers found in Seurat object\n")
    return(invisible(NULL))
  }
  
  unique_celltypes <- unique(seurat_obj@meta.data$temp_celltype_annotation)
  unique_celltypes <- unique_celltypes[!is.na(unique_celltypes) & unique_celltypes != "Unknown"]
  
  if (length(unique_celltypes) == 0) {
    cat("  ⚠ No annotated cell types found\n")
    return(invisible(NULL))
  }
  
  # Set identity
  Idents(seurat_obj) <- "temp_celltype_annotation"
  
  # Build DotPlot
  tryCatch({
    # Use Seurat::DotPlot to compute expression matrix, then rebuild with custom ggplot style
    p_base <- DotPlot(
      seurat_obj,
      features = valid_markers,
      dot.scale = 6
    )
    
    # Extract plot data (features.plot, id, pct.exp, avg.exp.scaled, etc.)
    data_plot <- p_base$data
    
    # Expression scaling range for fill gradient
    scale_min <- min(data_plot$avg.exp.scaled, na.rm = TRUE)
    scale_max <- max(data_plot$avg.exp.scaled, na.rm = TRUE)
    
    # Enforce factor levels for correct axis ordering
    # features.plot -> X axis (Markers)
    data_plot$features.plot <- factor(data_plot$features.plot, levels = valid_markers)
    # id -> Y axis (Cell Types)
    data_plot$id <- factor(data_plot$id, levels = ordered_celltypes)

    # Custom ggplot DotPlot:
    # - fill mapped to avg.exp.scaled with red gradient
    # - fixed black stroke (shape = 21: colored fill + black outline)
    p <- ggplot(
      data_plot,
      aes(x = features.plot, y = id)
    ) +
      ggrastr::rasterise(geom_point(
        shape = 21,
        colour = "black",
        stroke = 0.9,
        aes(size = pct.exp, fill = avg.exp.scaled)
      ), dpi = 300) +
      scale_size(range = c(0, 8)) +
      scale_fill_gradientn(
        colours = c("#fee0d2", "#fc9272", "#de2d26"),
        limits = c(scale_min, scale_max)
      ) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_line(color = adjustcolor("gray", alpha.f = 0.3))
      ) +
      labs(
        title = "Key Markers Validated by LLM",
        x = "Markers",
        y = "Cell Type"
      )
    
    # Calculate appropriate figure dimensions
    plot_width <- max(8, length(valid_markers) * 0.3 + 3)
    plot_height <- max(6, length(unique_celltypes) * 0.3 + 2)
    
    # Cap maximum dimensions
    plot_width <- min(plot_width, 20)
    plot_height <- min(plot_height, 15)
    
    # Save PDF
    showtext_auto(FALSE)
    ggsave(
      paste0(output_prefix, "_marker_dotplot.pdf"),
      p,
      width = plot_width,
      height = plot_height,
      device = cairo_pdf
    )
    showtext_auto(TRUE)
    
    # Save PNG
    ggsave(
      paste0(output_prefix, "_marker_dotplot.png"),
      p,
      width = plot_width,
      height = plot_height,
      dpi = 300,
      bg = "white"
    )
    
    cat(sprintf("  ✓ Saved: %s_marker_dotplot.pdf/png\n", basename(output_prefix)))
    cat(sprintf("    - Markers shown: %d (of %d requested)\n", 
                length(valid_markers), length(all_markers)))
    cat(sprintf("    - Cell types: %d\n", length(unique_celltypes)))
    
  }, error = function(e) {
    cat(sprintf("  ✗ DotPlot generation failed: %s\n", e$message))
  })
  
  # Remove temporary column
  seurat_obj@meta.data$temp_celltype_annotation <- NULL
  
  invisible(NULL)
}


cat("✓ Visualization functions loaded\n")
