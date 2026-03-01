# ============================================================================
# Data Preparation Functions
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})


#' Extract top-N candidates from ScType score matrix
#' 
#' @param score_matrix ScType score matrix (cell types x cells)
#' @param seurat_obj Seurat object
#' @param cluster_column Column name in meta.data containing cluster assignments
#' @param top_n Number of top candidates to extract per cluster
#' @return data.frame with top-N candidates for all clusters
extract_top_n_candidates <- function(
  score_matrix,
  seurat_obj,
  cluster_column = "seurat_clusters",
  top_n = 10
) {
  
  clusters <- unique(seurat_obj@meta.data[[cluster_column]])
  
  results_list <- list()
  
  for (cl in clusters) {
    cluster_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[cluster_column]] == cl, ])
    valid_cells <- intersect(cluster_cells, colnames(score_matrix))
    
    if (length(valid_cells) == 0) next
    
    es_subset <- score_matrix[, valid_cells, drop = FALSE]
    
    if (nrow(es_subset) == 1) {
      scores <- if (ncol(es_subset) == 1) {
        setNames(es_subset[1, 1], rownames(es_subset)[1])
      } else {
        setNames(sum(es_subset[1, ]), rownames(es_subset)[1])
      }
    } else {
      scores <- rowSums(es_subset)
    }
    
    scores_sorted <- sort(scores, decreasing = TRUE)
    top_scores <- head(scores_sorted, top_n)
    
    results_list[[as.character(cl)]] <- data.frame(
      cluster = cl,
      type = names(top_scores),
      scores = as.numeric(top_scores),
      ncells = sum(seurat_obj@meta.data[[cluster_column]] == cl),
      rank = seq_along(top_scores),
      stringsAsFactors = FALSE
    )
  }
  
  result_df <- do.call(rbind, results_list)
  rownames(result_df) <- NULL
  
  return(result_df)
}


#' Prepare Seurat object (normalize / variable features / scale / PCA / neighbors / optional clustering / optional UMAP)
#' 
#' Design goals:
#' 1) Build a single neighbor graph (fixed dims/k/assay) that can be reused across resolutions;
#' 2) Follow the canonical Seurat workflow: Normalize -> FindVariableFeatures -> Scale -> PCA -> FindNeighbors -> (optional) FindClusters -> (optional) UMAP.
#' 
#' @param seurat_obj Seurat object
#' @param force_reprocess If TRUE, forces full reprocessing from scratch
#' @param cluster_resolution Initial clustering resolution (only used when compute_clusters = TRUE)
#' @param assay Assay to process (default "RNA")
#' @param dims PCA dimensions for neighbor graph and UMAP, e.g. 1:30
#' @param k_param k parameter for KNN (passed to FindNeighbors k.param)
#' @param compute_clusters Whether to run FindClusters inside this function (default TRUE)
#' @param random_seed Random seed for reproducible clustering
#' @param run_umap Whether to compute UMAP embedding (default TRUE)
#' @return Processed Seurat object
prepare_seurat_object <- function(
  seurat_obj,
  force_reprocess = FALSE,
  cluster_resolution = 0.8,
  assay = "RNA",
  dims = 1:30,
  k_param = 20,
  compute_clusters = TRUE,
  random_seed = 1234,
  run_umap = TRUE
) {
  
  DefaultAssay(seurat_obj) <- assay
  
  has_pca <- "pca" %in% names(seurat_obj@reductions)
  has_graphs <- length(seurat_obj@graphs) > 0
  
  if (force_reprocess) {
    cat("  Processing Seurat object (forced rebuild)...\n")
    
    cat("    - Normalizing...\n")
    seurat_obj <- NormalizeData(seurat_obj, assay = assay, verbose = FALSE)
    
    cat("    - Finding variable features...\n")
    seurat_obj <- FindVariableFeatures(seurat_obj, assay = assay, verbose = FALSE)
    
    cat("    - Scaling variable features...\n")
    seurat_obj <- ScaleData(
      seurat_obj,
      assay = assay,
      features = VariableFeatures(seurat_obj),
      verbose = FALSE
    )
    
    cat("    - Running PCA...\n")
    seurat_obj <- RunPCA(seurat_obj, assay = assay, features = VariableFeatures(seurat_obj), verbose = FALSE)
    
    cat(sprintf("    - Building neighbors (dims=%s, k=%d)...\n", paste(range(dims), collapse = ":"), k_param))
    seurat_obj <- FindNeighbors(seurat_obj, dims = dims, k.param = k_param, verbose = FALSE)
    
    if (compute_clusters) {
      cat(sprintf("    - Clustering (resolution = %.3f)...\n", cluster_resolution))
      set.seed(random_seed)
      seurat_obj <- FindClusters(
        seurat_obj,
        resolution = cluster_resolution,
        random.seed = random_seed,
        verbose = FALSE
      )
    }
    
    if (run_umap) {
      cat(sprintf("    - Running UMAP (dims=%s)...\n", paste(range(dims), collapse = ":")))
      seurat_obj <- RunUMAP(seurat_obj, dims = dims, verbose = FALSE)
    }
    
    cat("  ✓ Processing complete\n")
    return(seurat_obj)
  }
  
  # Non-forced: fill in only missing steps without rebuilding the existing neighbor graph,
  # allowing the same graph to be reused across multiple clustering resolutions.
  if (!has_pca) {
    cat("  Processing Seurat object (filling missing steps)...\n")
    
    cat("    - Normalizing...\n")
    seurat_obj <- NormalizeData(seurat_obj, assay = assay, verbose = FALSE)
    
    cat("    - Finding variable features...\n")
    seurat_obj <- FindVariableFeatures(seurat_obj, assay = assay, verbose = FALSE)
    
    cat("    - Scaling variable features...\n")
    seurat_obj <- ScaleData(
      seurat_obj,
      assay = assay,
      features = VariableFeatures(seurat_obj),
      verbose = FALSE
    )
    
    cat("    - Running PCA...\n")
    seurat_obj <- RunPCA(seurat_obj, assay = assay, features = VariableFeatures(seurat_obj), verbose = FALSE)
  }
  
  if (!has_graphs) {
    cat(sprintf("    - Building neighbors (dims=%s, k=%d)...\n", paste(range(dims), collapse = ":"), k_param))
    seurat_obj <- FindNeighbors(seurat_obj, dims = dims, k.param = k_param, verbose = FALSE)
  }
  
  if (compute_clusters && !"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    cat(sprintf("    - Clustering (resolution = %.3f)...\n", cluster_resolution))
    set.seed(random_seed)
    seurat_obj <- FindClusters(
      seurat_obj,
      resolution = cluster_resolution,
      random.seed = random_seed,
      verbose = FALSE
    )
  }
  
  if (run_umap && !"umap" %in% names(seurat_obj@reductions)) {
    cat(sprintf("    - Running UMAP (dims=%s)...\n", paste(range(dims), collapse = ":")))
    seurat_obj <- RunUMAP(seurat_obj, dims = dims, verbose = FALSE)
  }
  
  cat("  ✓ Seurat object ready\n")
  return(seurat_obj)
}


#' Extract scaled data (compatible with Seurat v4 and v5)
#' 
#' @param seurat_obj Seurat object
#' @return Scaled data matrix
extract_scaled_data <- function(seurat_obj) {
  
  seurat_v5 <- isFALSE('counts' %in% names(attributes(seurat_obj[["RNA"]])))
  
  if (seurat_v5) {
    scaled_data <- as.matrix(seurat_obj[["RNA"]]$scale.data)
  } else {
    scaled_data <- as.matrix(seurat_obj[["RNA"]]@scale.data)
  }
  
  return(scaled_data)
}


#' Generate final annotation summary
#' 
#' @param llm_results LLM decision results
#' @param sctype_annotation Original ScType annotation
#' @param hypergeom_results Hypergeometric test results
#' @return data.frame containing the final annotation summary
generate_final_annotation_summary <- function(
  llm_results,
  sctype_annotation,
  hypergeom_results
) {
  
  candidate_stats <- hypergeom_results %>%
    group_by(cluster) %>%
    summarise(
      n_candidates_total = n(),
      n_candidates_significant = sum(is_significant),
      top_sctype_score = max(sctype_score),
      .groups = "drop"
    )
  
  final_summary <- llm_results %>%
    left_join(candidate_stats, by = "cluster") %>%
    left_join(
      sctype_annotation %>% select(cluster, sctype_celltype = type),
      by = "cluster"
    ) %>%
    mutate(
      changed_from_sctype = selected_celltype != sctype_celltype
    ) %>%
    select(
      cluster,
      llm_celltype = selected_celltype,
      confidence,
      sctype_celltype,
      changed_from_sctype,
      n_candidates_total,
      n_candidates_significant,
      key_markers_validated,
      reasoning
    )
  
  return(final_summary)
}

cat("✓ Data preparation functions loaded\n")
