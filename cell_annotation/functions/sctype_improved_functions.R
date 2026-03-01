#!/usr/bin/env Rscript
# ============================================================================
# ScType Improved Algorithm - Functions for MyData Database
# ============================================================================
# Modified version of ScType to work with custom mydata database format
# Features:
#   - Species filtering
#   - Tissue filtering (top_level + top_level_mixed_name)
#   - Quality control filtering (gene_qc_pass)
#   - Support for positive/negative marker polarity
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(HGNChelper)
  library(scales)
})

# ============================================================================
# Function: prepare_gene_sets_from_mydata
# ============================================================================
# Prepare gene sets from mydata database with filtering options
#
# Parameters:
#   mydata: data.table/data.frame with marker database
#   species: Species filter (e.g., "Human", "Mouse", NULL for all)
#   tissue: Tissue filter using top_level (NULL for all)
#   tissue_mixed_name: Optional specific mixed tissue name filter
#   use_standardized_cell_name: Use cell_name_standardized (TRUE) or cell_name_cl (FALSE)
#   min_markers_per_cell: Minimum number of markers required per cell type
#   use_qc_filter: Only use markers with gene_qc_pass == TRUE
#   disease_filter: Character scalar/vector for disease_type filtering. Examples:
#                   "Normal" (default) keeps only normal entries;
#                   c("breast cancer", "lung cancer") keeps only these diseases;
#                   NULL disables disease filtering.
#   disease_column: Column name to use for disease filtering (default: "disease_type")
#                   Can be set to "disease_type_do" for DO-based filtering
#   use_normal_cells_only: Backward-compat; if TRUE and disease_filter is NULL, use "Normal"
#   verbose: Print filtering statistics
#
# Returns:
#   List with:
#     - gs_positive: List of positive marker gene sets per cell type
#     - gs_negative: List of negative marker gene sets per cell type
#     - metadata: Filtering statistics
#
prepare_gene_sets_from_mydata <- function(
    mydata,
    species = "Human",
    tissue = NULL,
    tissue_mixed_name = NULL,
    use_standardized_cell_name = TRUE,
    min_markers_per_cell = 1,
    use_qc_filter = TRUE,
    disease_filter = NULL,
    disease_column = "disease_type",
    use_normal_cells_only = TRUE,
    verbose = TRUE
) {
  
  if (verbose) {
    cat("\n")
    cat("========================================================================\n")
    cat("  Preparing Gene Sets from MyData Database\n")
    cat("========================================================================\n\n")
    cat(sprintf("Original database size: %d entries\n", nrow(mydata)))
  }
  
  # Convert to data.table if needed
  if (!inherits(mydata, "data.table")) {
    mydata <- as.data.table(mydata)
  }
  
  # Make a copy to avoid modifying original
  db_filtered <- copy(mydata)
  
  # -------------------------------------------------------------------------
  # Step 1: Filter by species
  # -------------------------------------------------------------------------
  if (!is.null(species)) {
    # Avoid data.table NSE capturing the 'species' column instead of the function arg
    species_filter_values <- as.character(species)
    db_filtered <- db_filtered[as.character(db_filtered$species) %in% species_filter_values, ]
    if (verbose) {
      cat(sprintf("After species filter (%s): %d entries\n",
                  paste(species_filter_values, collapse = ", "), nrow(db_filtered)))
    }
  }
  
  if (nrow(db_filtered) == 0) {
    stop("No entries remain after species filtering!")
  }
  
  # -------------------------------------------------------------------------
  # Step 2: Filter by tissue
  # -------------------------------------------------------------------------
  if (!is.null(tissue)) {
    # Handle mixed tissue separately if specified
    if (!is.null(tissue_mixed_name)) {
      # For mixed_tissue, also filter by mixed_name
      db_filtered <- db_filtered[
        (db_filtered$top_level %in% tissue & db_filtered$top_level != "mixed_tissue") |
        (db_filtered$top_level == "mixed_tissue" & db_filtered$top_level_mixed_name %in% tissue_mixed_name),
        
      ]
      if (verbose) {
        cat(sprintf("After tissue filter (top_level: %s, mixed_name: %s): %d entries\n",
                    paste(tissue, collapse = ", "),
                    paste(tissue_mixed_name, collapse = ", "),
                    nrow(db_filtered)))
      }
    } else {
      # Normal tissue filtering, exclude mixed_tissue unless explicitly requested
      if ("mixed_tissue" %in% tissue) {
        db_filtered <- db_filtered[db_filtered$top_level %in% tissue, ]
      } else {
        db_filtered <- db_filtered[db_filtered$top_level %in% tissue & db_filtered$top_level != "mixed_tissue", ]
      }
      if (verbose) {
        cat(sprintf("After tissue filter (%s): %d entries\n",
                    paste(tissue, collapse = ", "), nrow(db_filtered)))
      }
    }
  }
  
  if (nrow(db_filtered) == 0) {
    stop("No entries remain after tissue filtering!")
  }
  
  # -------------------------------------------------------------------------
  # Step 3: Filter by disease using specified disease_column (preferred over use_normal_cells_only)
  # -------------------------------------------------------------------------
  # Backward-compat: if disease_filter is NULL and use_normal_cells_only is TRUE, set to "Normal"
  if (is.null(disease_filter) && isTRUE(use_normal_cells_only)) {
    disease_filter <- "Normal"
  }

  if (!is.null(disease_filter)) {
    if (!(disease_column %in% colnames(db_filtered))) {
      stop(sprintf("disease_filter is set but '%s' column not found in mydata", disease_column))
    }

    # Ensure vector
    disease_filter_vec <- as.character(disease_filter)

    # Exact match filtering using the specified disease column
    db_filtered <- db_filtered[get(disease_column) %in% disease_filter_vec]

    if (verbose) {
      cat(sprintf("After %s filter (%s): %d entries\n",
                  disease_column,
                  paste(disease_filter_vec, collapse = ", "), nrow(db_filtered)))
    }
  }
  
  # -------------------------------------------------------------------------
  # Step 4: Filter by gene QC status
  # -------------------------------------------------------------------------
  if (use_qc_filter) {
    # Convert gene_qc_pass to logical if it's character
    if (is.character(db_filtered$gene_qc_pass)) {
      db_filtered$gene_qc_pass <- db_filtered$gene_qc_pass == "TRUE"
    }
    db_filtered <- db_filtered[gene_qc_pass == TRUE | is.na(gene_qc_pass)]
    if (verbose) {
      cat(sprintf("After gene QC filter: %d entries\n", nrow(db_filtered)))
    }
  }
  
  # -------------------------------------------------------------------------
  # Step 5: Filter entries with valid gene symbols
  # -------------------------------------------------------------------------
  db_filtered <- db_filtered[!is.na(gene_symbol) & gene_symbol != ""]
  if (verbose) {
    cat(sprintf("After removing entries without gene symbols: %d entries\n", nrow(db_filtered)))
  }
  
  if (nrow(db_filtered) == 0) {
    stop("No entries remain after all filtering steps!")
  }
  
  # -------------------------------------------------------------------------
  # Step 6: Choose cell name column
  # -------------------------------------------------------------------------
  if (use_standardized_cell_name) {
    cell_name_col <- "cell_name_standardized"
  } else {
    cell_name_col <- "cell_name_cl"
  }
  
  # Remove entries with missing cell names
  db_filtered <- db_filtered[!is.na(get(cell_name_col)) & get(cell_name_col) != ""]
  if (verbose) {
    cat(sprintf("Using cell name column: %s\n", cell_name_col))
    cat(sprintf("After removing entries without cell names: %d entries\n", nrow(db_filtered)))
  }
  
  # -------------------------------------------------------------------------
  # Step 7: Separate positive and negative markers
  # -------------------------------------------------------------------------
  db_positive <- db_filtered[marker_polarity == "positive" | marker_polarity == "positive"]
  db_negative <- db_filtered[marker_polarity == "negative" | marker_polarity == "negative"]
  
  if (verbose) {
    cat(sprintf("Positive markers: %d entries\n", nrow(db_positive)))
    cat(sprintf("Negative markers: %d entries\n", nrow(db_negative)))
  }
  
  # -------------------------------------------------------------------------
  # Step 8: Aggregate markers by cell type
  # -------------------------------------------------------------------------
  
  # Process positive markers
  positive_markers_by_cell <- db_positive %>%
    group_by(cell_type = .data[[cell_name_col]]) %>%
    summarise(
      markers = list(unique(toupper(gene_symbol))),
      n_markers = n_distinct(gene_symbol),
      .groups = "drop"
    ) %>%
    filter(n_markers >= min_markers_per_cell)
  
  # Process negative markers
  negative_markers_by_cell <- db_negative %>%
    group_by(cell_type = .data[[cell_name_col]]) %>%
    summarise(
      markers = list(unique(toupper(gene_symbol))),
      n_markers = n_distinct(gene_symbol),
      .groups = "drop"
    )
  
  if (verbose) {
    cat(sprintf("\nCell types with >= %d markers: %d\n", 
                min_markers_per_cell, nrow(positive_markers_by_cell)))
    cat(sprintf("Cell types with negative markers: %d\n", nrow(negative_markers_by_cell)))
  }
  
  # -------------------------------------------------------------------------
  # Step 9: Create gene sets lists
  # -------------------------------------------------------------------------
  
  # Positive markers
  gs_positive <- setNames(
    positive_markers_by_cell$markers,
    positive_markers_by_cell$cell_type
  )
  
  # Negative markers (fill with empty vector for cell types without negative markers)
  gs_negative <- setNames(
    lapply(positive_markers_by_cell$cell_type, function(ct) {
      neg_idx <- which(negative_markers_by_cell$cell_type == ct)
      if (length(neg_idx) > 0) {
        negative_markers_by_cell$markers[[neg_idx]]
      } else {
        character(0)
      }
    }),
    positive_markers_by_cell$cell_type
  )
  
  # -------------------------------------------------------------------------
  # Step 10: Skip gene symbol correction
  # -------------------------------------------------------------------------
  # Keep marker gene symbols as provided; only ensure alignment and non-empty sets.
  
  # Remove cell types with no markers
  valid_idx <- sapply(gs_positive, length) > 0
  gs_positive <- gs_positive[valid_idx]
  gs_negative <- gs_negative[valid_idx]
  
  if (verbose) {
    cat(sprintf("Final cell types: %d\n", length(gs_positive)))
    
    # Print summary statistics
    marker_counts <- sapply(gs_positive, length)
    cat(sprintf("\nMarker statistics (positive):\n"))
    cat(sprintf("  Mean markers per cell type: %.1f\n", mean(marker_counts)))
    cat(sprintf("  Median markers per cell type: %.0f\n", median(marker_counts)))
    cat(sprintf("  Range: %d - %d\n", min(marker_counts), max(marker_counts)))
    
    cat(sprintf("\nTop 10 cell types by marker count:\n"))
    top_cells <- sort(marker_counts, decreasing = TRUE)[1:min(10, length(marker_counts))]
    for (i in seq_along(top_cells)) {
      cat(sprintf("  %2d. %-50s: %3d markers\n", i, names(top_cells)[i], top_cells[i]))
    }
    
    cat("\n========================================================================\n\n")
  }
  
  # -------------------------------------------------------------------------
  # Step 11: Create metadata
  # -------------------------------------------------------------------------
  metadata <- list(
    original_entries = nrow(mydata),
    filtered_entries = nrow(db_filtered),
    species_filter = species,
    tissue_filter = tissue,
    tissue_mixed_name_filter = tissue_mixed_name,
    disease_filter = disease_filter,
    cell_name_column = cell_name_col,
    n_cell_types = length(gs_positive),
    n_positive_markers = sum(sapply(gs_positive, length)),
    n_negative_markers = sum(sapply(gs_negative, length)),
    use_qc_filter = use_qc_filter,
    use_normal_cells_only = use_normal_cells_only
  )
  
  return(list(
    gs_positive = gs_positive,
    gs_negative = gs_negative,
    metadata = metadata
  ))
}


# ============================================================================
# Function: sctype_score
# ============================================================================
# Calculate ScType scores for cell type assignment
# (Adapted from original ScType algorithm)
#
# Parameters:
#   scRNAseqData: Gene expression matrix (genes x cells)
#   scaled: Whether the data is already scaled
#   gs: List of positive marker genes per cell type
#   gs2: List of negative marker genes per cell type (optional)
#   gene_names_to_uppercase: Convert gene names to uppercase
#   verbose: Print progress messages
#
sctype_score <- function(scRNAseqData, scaled = TRUE, gs, gs2 = NULL, 
                         gene_names_to_uppercase = TRUE, verbose = TRUE) {
  
  # Check input matrix
  if (!is.matrix(scRNAseqData)) {
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if (sum(dim(scRNAseqData)) == 0) {
      warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # Handle NULL gs2
  if (is.null(gs2)) {
    gs2 <- lapply(gs, function(x) character(0))
    names(gs2) <- names(gs)
  }
  
  if (verbose) {
    cat(sprintf("Calculating ScType scores for %d cell types across %d cells...\n",
                length(gs), ncol(scRNAseqData)))
  }
  
  # Marker sensitivity (genes appearing in fewer cell types get higher weight)
  marker_stat <- sort(table(unlist(gs)), decreasing = TRUE)
  marker_sensitivity <- data.frame(
    score_marker_sensitivity = scales::rescale(
      as.numeric(marker_stat), 
      to = c(0, 1), 
      from = c(length(gs), 1)
    ),
    gene_ = names(marker_stat), 
    stringsAsFactors = FALSE
  )
  
  # Convert gene names to uppercase
  if (gene_names_to_uppercase) {
    rownames(scRNAseqData) <- toupper(rownames(scRNAseqData))
  }
  
  # Subselect genes only found in data
  names_gs_cp <- names(gs)
  names_gs_2_cp <- names(gs2)
  
  gs <- lapply(1:length(gs), function(d_) { 
    GeneIndToKeep <- rownames(scRNAseqData) %in% as.character(gs[[d_]])
    rownames(scRNAseqData)[GeneIndToKeep]
  })
  
  gs2 <- lapply(1:length(gs2), function(d_) { 
    GeneIndToKeep <- rownames(scRNAseqData) %in% as.character(gs2[[d_]])
    rownames(scRNAseqData)[GeneIndToKeep]
  })
  
  names(gs) <- names_gs_cp
  names(gs2) <- names_gs_2_cp
  
  cell_markers_genes_score <- marker_sensitivity[
    marker_sensitivity$gene_ %in% unique(unlist(gs)), 
  ]
  
  # Z-scale if not already scaled
  if (!scaled) {
    Z <- t(scale(t(scRNAseqData)))
  } else {
    Z <- scRNAseqData
  }
  
  # Multiply by marker sensitivity
  for (jj in 1:nrow(cell_markers_genes_score)) {
    gene_name <- cell_markers_genes_score[jj, "gene_"]
    if (gene_name %in% rownames(Z)) {
      Z[gene_name, ] <- Z[gene_name, ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
    }
  }
  
  # Subselect only marker genes
  all_markers <- unique(c(unlist(gs), unlist(gs2)))
  Z <- Z[rownames(Z) %in% all_markers, , drop = FALSE]
  
  if (verbose) {
    cat(sprintf("Using %d marker genes found in the data\n", nrow(Z)))
  }
  
  # ============================================================================
  # Calculate enrichment scores
  # ============================================================================
  
  es <- do.call("rbind", lapply(names(gs), function(gss_) { 
    sapply(1:ncol(Z), function(j) {
      gs_z <- Z[gs[[gss_]], j]
      gz_2 <- Z[gs2[[gss_]], j] * -1
      
      sum_t1 <- if (length(gs_z) > 0) {
        sum(gs_z) / sqrt(length(gs_z))
      } else {
        0
      }
      
      sum_t2 <- if (length(gz_2) > 0) {
        sum(gz_2) / sqrt(length(gz_2))
      } else {
        0
      }
      
      if (is.na(sum_t2)) sum_t2 <- 0
      
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) <- list(names(gs), colnames(Z))
  
  # Remove NA rows
  es.max <- es[!apply(is.na(es) | es == "", 1, all), , drop = FALSE]
  
  if (verbose) {
    cat(sprintf("Score matrix: %d cell types x %d cells\n", 
                nrow(es.max), ncol(es.max)))
  }
  
  return(es.max)
}


# ============================================================================
# Function: assign_cell_types_by_cluster
# ============================================================================
# Assign cell types to Seurat clusters based on ScType scores
#
# Parameters:
#   seurat_obj: Seurat object with clustering
#   score_matrix: ScType score matrix (cell types x cells)
#   cluster_column: Column name in meta.data containing cluster assignments
#   min_score_ratio: Minimum ratio of score to ncells for confident assignment
#   return_all_scores: Return top N scores per cluster (for inspection)
#
assign_cell_types_by_cluster <- function(
    seurat_obj,
    score_matrix,
    cluster_column = "seurat_clusters",
    min_score_ratio = 0.25,
    return_all_scores = TRUE,
    top_n_per_cluster = 10,
    verbose = TRUE
) {
  
  if (verbose) {
    cat("\nAssigning cell types to clusters...\n")
  }
  
  clusters <- unique(seurat_obj@meta.data[[cluster_column]])
  
  # Calculate scores for each cluster
  cL_results <- do.call("rbind", lapply(clusters, function(cl) {
    cluster_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[cluster_column]] == cl, ])
    valid_cells <- intersect(cluster_cells, colnames(score_matrix))
    
    if (length(valid_cells) == 0) return(NULL)
    
    es.max.cl <- tryCatch({
      es_subset <- score_matrix[, valid_cells, drop = FALSE]
      if (!is.matrix(es_subset) || nrow(es_subset) == 0) return(NULL)
      
      if (nrow(es_subset) == 1) {
        scores <- if (ncol(es_subset) == 1) {
          setNames(es_subset[1, 1], rownames(es_subset)[1])
        } else {
          setNames(sum(es_subset[1, ]), rownames(es_subset)[1])
        }
        sort(scores, decreasing = TRUE)
      } else {
        sort(rowSums(es_subset), decreasing = TRUE)
      }
    }, error = function(e) NULL)
    
    if (is.null(es.max.cl) || length(es.max.cl) == 0) return(NULL)
    
    head(data.frame(
      cluster = cl,
      type = names(es.max.cl),
      scores = as.numeric(es.max.cl),
      ncells = sum(seurat_obj@meta.data[[cluster_column]] == cl),
      stringsAsFactors = FALSE
    ), top_n_per_cluster)
  }))
  
  # Get top cell type per cluster
  sctype_scores <- cL_results %>%
    group_by(cluster) %>%
    slice_max(order_by = scores, n = 1) %>%
    ungroup()
  
  # Mark low-confidence predictions as "Unknown"
  sctype_scores$type[sctype_scores$scores < sctype_scores$ncells * min_score_ratio] <- "Unknown"
  
  if (verbose) {
    cat(sprintf("Annotated %d clusters\n", nrow(sctype_scores)))
    cat(sprintf("Confident assignments: %d\n", sum(sctype_scores$type != "Unknown")))
    cat(sprintf("Unknown assignments: %d\n", sum(sctype_scores$type == "Unknown")))
  }
  
  # Create result list
  result <- list(
    cluster_annotation = sctype_scores
  )
  
  if (return_all_scores) {
    result$all_scores <- cL_results
  }
  
  return(result)
}


# ============================================================================
# Function: add_sctype_annotation_to_seurat
# ============================================================================
# Add ScType annotations to Seurat object metadata
#
add_sctype_annotation_to_seurat <- function(
    seurat_obj,
    cluster_annotation,
    cluster_column = "seurat_clusters",
    annotation_column = "sctype_annotation"
) {
  
  # Initialize annotation column
  seurat_obj@meta.data[[annotation_column]] <- ""
  
  # Assign annotations
  for (j in unique(cluster_annotation$cluster)) {
    cl_type <- cluster_annotation[cluster_annotation$cluster == j, ]
    seurat_obj@meta.data[[annotation_column]][
      seurat_obj@meta.data[[cluster_column]] == j
    ] <- as.character(cl_type$type[1])
  }
  
  return(seurat_obj)
}

cat("\n✓ ScType improved functions loaded successfully\n\n")

