#!/usr/bin/env Rscript
# ============================================================================
# scMarkerAgent - Cell Annotation Agent
# ============================================================================
# Automated cell type annotation using scMarkerAgent database and LLM validation
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(httr)
  library(jsonlite)
  library(ggplot2)
  library(patchwork)
  library(showtext)
  library(ggrastr)
})

# ----------------------------------------------------------------------------
# UTF-8 / font configuration for ggplot2
# Ensures that UTF-8 text (e.g., Greek letters like α, β) is rendered correctly
# in all generated plots.
# ----------------------------------------------------------------------------
font_add("Arial", "/usr/share/fonts/msttcore/arial.ttf")
showtext_auto()
showtext_opts(dpi = 300)

# ============================================================================
# Parse Command Line Arguments
# ============================================================================

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input Seurat object path (RDS file) [REQUIRED]", metavar = "FILE"),
  
  make_option("--n_variable_features", type = "integer", default = 2000,
              help = "Number of variable features [default: %default]"),
  
  make_option("--dims", type = "integer", default = 30,
              help = "PCA dimensions (1:N) [default: 1:%default]"),
  
  make_option("--cluster_resolutions", type = "character", default = "0.5",
              help = "Clustering resolutions, comma-separated (max 3) [default: %default]"),
  
  make_option("--random_seed", type = "integer", default = 1234,
              help = "Random seed for reproducibility [default: %default]"),
  
  make_option("--species", type = "character", default = "Human",
              help = "Species filter: Human, Mouse, or Rat [default: %default]"),
  
  make_option("--tissue", type = "character", default = NULL,
              help = "Tissue types, comma-separated [REQUIRED]", metavar = "STRING"),
  
  make_option("--min_markers_per_cell", type = "integer", default = 2,
              help = "Minimum markers per cell type [default: %default]"),
  
  make_option("--disease_type", type = "character", default = "Normal",
              help = "Disease types, comma-separated [default: %default]"),
  
  make_option(c("-o", "--output_dir"), type = "character", default = "./results",
              help = "Output directory [default: %default]"),
  
  make_option("--help_extended", action = "store_true", default = FALSE,
              help = "Show extended help with examples")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "\nscMarkerAgent - Cell Annotation Agent",
  epilogue = "For extended help with examples, use --help_extended"
)

opt <- parse_args(opt_parser)

# Extended help
if (opt$help_extended) {
  cat("\n")
  cat("============================================================================\n")
  cat("  scMarkerAgent - Extended Help\n")
  cat("============================================================================\n\n")
  cat("DESCRIPTION:\n")
  cat("  This pipeline performs automated cell type annotation using:\n")
  cat("  1. Standard Seurat preprocessing workflow\n")
  cat("  2. MyData database for cell type marker genes\n")
  cat("  3. ScType scoring algorithm\n")
  cat("  4. LLM validation for final annotation\n\n")
  
  cat("EXAMPLES:\n\n")
  cat("  # Basic usage with brain tissue:\n")
  cat("  Rscript annotate_cells.R \\\n")
  cat("    --input ./data/example_pbmc_input.rds \\\n")
  cat("    --tissue brain \\\n")
  cat("    --output_dir ./results/brain_annotation\n\n")
  
  cat("  # Multiple tissues and resolutions:\n")
  cat("  Rscript annotate_cells.R \\\n")
  cat("    --input ./data/sample.rds \\\n")
  cat("    --tissue \"brain,blood\" \\\n")
  cat("    --cluster_resolutions \"0.3,0.5,1.0\" \\\n")
  cat("    --n_variable_features 3000 \\\n")
  cat("    --dims 40\n\n")
  
  cat("  # Disease-specific annotation:\n")
  cat("  Rscript annotate_cells.R \\\n")
  cat("    --input ./data/cancer_sample.rds \\\n")
  cat("    --tissue lung \\\n")
  cat("    --disease_type \"Normal,lung cancer\"\n\n")
  
  cat("OUTPUT STRUCTURE:\n")
  cat("  results/\n")
  cat("  ├── annotated_seurat.RDS           # Final annotated Seurat object\n")
  cat("  ├── summary.txt                    # Run summary\n")
  cat("  └── resolution_X.X/                # Results per resolution\n")
  cat("      ├── cluster_markers_all.csv\n")
  cat("      ├── cluster_markers_significant.csv\n")
  cat("      ├── llm_decisions.csv\n")
  cat("      ├── llm_prompts_log.csv\n")
  cat("      ├── umap_cluster_annotation.pdf\n")
  cat("      └── umap_marker_dotplot.pdf\n\n")
  
  cat("FIXED PARAMETERS:\n")
  cat("  - logfc_threshold: 1\n")
  cat("  - top_n_candidates: 10\n")
  cat("  - assay: RNA\n")
  cat("  - k_param: 20\n\n")
  
  cat("============================================================================\n\n")
  quit(save = "no", status = 0)
}

# Validate required arguments
if (is.null(opt$input)) {
  cat("ERROR: --input is required\n\n")
  print_help(opt_parser)
  quit(save = "no", status = 1)
}

if (is.null(opt$tissue)) {
  cat("ERROR: --tissue is required\n\n")
  print_help(opt_parser)
  quit(save = "no", status = 1)
}

# Check input file exists
if (!file.exists(opt$input)) {
  cat(sprintf("ERROR: Input file not found: %s\n", opt$input))
  quit(save = "no", status = 1)
}

# Check file size (<= 5GB)
file_size_gb <- file.info(opt$input)$size / (1024^3)
if (file_size_gb > 5) {
  cat(sprintf("ERROR: Input file exceeds 5GB limit (%.2f GB)\n", file_size_gb))
  quit(save = "no", status = 1)
}

# Parse cluster resolutions
resolution_strs <- trimws(strsplit(opt$cluster_resolutions, ",")[[1]])
resolutions <- as.numeric(resolution_strs)
if (length(resolutions) > 3) {
  cat("ERROR: Maximum 3 cluster resolutions allowed\n")
  quit(save = "no", status = 1)
}

# Parse tissue list
tissues <- trimws(strsplit(opt$tissue, ",")[[1]])

# Parse disease types
disease_types <- trimws(strsplit(opt$disease_type, ",")[[1]])

# Create output directory
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# ============================================================================
# Configuration
# ============================================================================

# Set LLM API configuration
Sys.setenv(PRIMARY_API_KEY = "YOUR_API_KEY_HERE")
Sys.setenv(PRIMARY_API_URL = "https://api.openai.com/v1/chat/completions")
Sys.setenv(PRIMARY_MODEL = "gpt-5-2025-08-07")

Sys.setenv(SECONDARY_API_KEY = "YOUR_API_KEY_HERE")
Sys.setenv(SECONDARY_API_URL = "https://api.openai.com/v1/chat/completions")
Sys.setenv(SECONDARY_MODEL = "gpt-5-2025-08-07")

# Fixed parameters
LOGFC_THRESHOLD <- 1
TOP_N_CANDIDATES <- 10
ASSAY <- "RNA"
K_PARAM <- 20

# Load functions
get_script_dir <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    return(dirname(normalizePath(sub(needle, "", cmdArgs[match]))))
  } else {
    return(getwd())
  }
}

script_dir <- get_script_dir()
source(file.path(script_dir, "functions", "sctype_improved_functions.R"))
source(file.path(script_dir, "functions", "data_preparation.R"))
source(file.path(script_dir, "functions", "llm_api.R"))
source(file.path(script_dir, "functions", "visualization.R"))

# ============================================================================
# Print Configuration
# ============================================================================

cat("\n")
cat("========================================================================\n")
cat("  scMarkerAgent - Cell Annotation Agent\n")
cat("========================================================================\n\n")

cat("Configuration:\n")
cat(sprintf("  Input file: %s (%.2f GB)\n", opt$input, file_size_gb))
cat(sprintf("  Output directory: %s\n", opt$output_dir))
cat(sprintf("  Species: %s\n", opt$species))
cat(sprintf("  Tissues: %s\n", paste(tissues, collapse = ", ")))
cat(sprintf("  Disease types: %s\n", paste(disease_types, collapse = ", ")))
cat(sprintf("  Variable features: %d\n", opt$n_variable_features))
cat(sprintf("  PCA dimensions: 1:%d\n", opt$dims))
cat(sprintf("  Cluster resolutions: %s\n", paste(resolutions, collapse = ", ")))
cat(sprintf("  Random seed: %d\n", opt$random_seed))
cat("\nFixed parameters:\n")
cat(sprintf("  LogFC threshold: %d\n", LOGFC_THRESHOLD))
cat(sprintf("  Top N candidates: %d\n", TOP_N_CANDIDATES))
cat(sprintf("  K-param (KNN): %d\n", K_PARAM))
cat("\n")

start_time <- Sys.time()

# ============================================================================
# Step 1: Load Input Data
# ============================================================================

cat("========================================================================\n")
cat("  Step 1: Loading Input Data\n")
cat("========================================================================\n\n")

cat("Loading Seurat object...\n")
seurat_obj <- readRDS(opt$input)
cat(sprintf("  ✓ Loaded: %d cells, %d features\n", ncol(seurat_obj), nrow(seurat_obj)))

cat("\nLoading MyData database...\n")
mydata_path <- file.path(dirname(script_dir), "data", "mydata.RDS")
if (!file.exists(mydata_path)) {
  cat(sprintf("ERROR: MyData database not found: %s\n", mydata_path))
  quit(save = "no", status = 1)
}
mydata <- readRDS(mydata_path)
cat(sprintf("  ✓ Loaded: %d entries\n\n", nrow(mydata)))

# ============================================================================
# Step 2: Preprocessing
# ============================================================================

cat("========================================================================\n")
cat("  Step 2: Preprocessing\n")
cat("========================================================================\n\n")

DefaultAssay(seurat_obj) <- ASSAY

cat("Step 2.1: Normalizing data...\n")
seurat_obj <- NormalizeData(seurat_obj, assay = ASSAY, verbose = FALSE)

cat("Step 2.2: Finding variable features...\n")
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  assay = ASSAY,
  selection.method = "vst",
  nfeatures = opt$n_variable_features,
  verbose = FALSE
)

cat("Step 2.3: Scaling data...\n")
seurat_obj <- ScaleData(
  seurat_obj,
  assay = ASSAY,
  features = VariableFeatures(seurat_obj),
  verbose = FALSE
)

cat("Step 2.4: Running PCA...\n")
set.seed(opt$random_seed)
seurat_obj <- RunPCA(
  seurat_obj,
  assay = ASSAY,
  features = VariableFeatures(seurat_obj),
  seed.use = opt$random_seed,
  verbose = FALSE
)

cat("Step 2.5: Building neighborhood graph...\n")
seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = 1:opt$dims,
  k.param = K_PARAM,
  verbose = FALSE
)

cat("Step 2.6: Running UMAP...\n")
set.seed(opt$random_seed)
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:opt$dims,
  seed.use = opt$random_seed,
  verbose = FALSE
)

cat("\nStep 2.7: Clustering at multiple resolutions...\n")
for (res in resolutions) {
  cat(sprintf("  - Resolution %.2f...\n", res))
  set.seed(opt$random_seed)
  seurat_obj <- FindClusters(
    seurat_obj,
    resolution = res,
    random.seed = opt$random_seed,
    verbose = FALSE
  )
  cluster_col <- paste0("RNA_snn_res.", res)
  n_clusters <- length(unique(seurat_obj@meta.data[[cluster_col]]))
  cat(sprintf("    ✓ %d clusters\n", n_clusters))
}

cat("\n  ✓ Preprocessing completed\n\n")

# ============================================================================
# Step 3: Prepare Gene Sets from MyData
# ============================================================================

cat("========================================================================\n")
cat("  Step 3: Preparing Gene Sets from MyData\n")
cat("========================================================================\n\n")

gs_list <- prepare_gene_sets_from_mydata(
  mydata = mydata,
  species = opt$species,
  tissue = tissues,
  use_standardized_cell_name = TRUE,
  min_markers_per_cell = opt$min_markers_per_cell,
  use_qc_filter = FALSE,
  disease_filter = disease_types,
  disease_column = "disease_type_do",
  use_normal_cells_only = FALSE,
  verbose = TRUE
)

cat(sprintf("  ✓ Gene sets prepared: %d cell types\n\n", length(gs_list$gs_positive)))

# ============================================================================
# Step 4: Calculate ScType Scores
# ============================================================================

cat("========================================================================\n")
cat("  Step 4: Calculating ScType Scores\n")
cat("========================================================================\n\n")

scaled_data <- extract_scaled_data(seurat_obj)
cat(sprintf("Scaled data: %d genes x %d cells\n\n", nrow(scaled_data), ncol(scaled_data)))

es_max <- sctype_score(
  scRNAseqData = scaled_data,
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = gs_list$gs_negative,
  verbose = TRUE
)

cat(sprintf("\n  ✓ Score matrix: %d cell types x %d cells\n\n", nrow(es_max), ncol(es_max)))

# ============================================================================
# Step 5: Process Each Resolution with LLM Validation
# ============================================================================

cat("========================================================================\n")
cat("  Step 5: Processing Each Resolution with LLM Validation\n")
cat("========================================================================\n\n")

for (i in seq_along(resolutions)) {
  res <- resolutions[i]
  res_str <- resolution_strs[i]
  cat("----------------------------------------------------------------------\n")
  cat(sprintf("  Resolution: %s\n", res_str))
  cat("----------------------------------------------------------------------\n\n")
  
  cluster_col <- paste0("RNA_snn_res.", res)
  res_dir <- file.path(opt$output_dir, paste0("resolution_", res_str))
  if (!dir.exists(res_dir)) {
    dir.create(res_dir, recursive = TRUE)
  }
  
  # 5.1 Get ScType annotation and top candidates
  cat("Step 5.1: ScType annotation...\n")
  sctype_results <- assign_cell_types_by_cluster(
    seurat_obj = seurat_obj,
    score_matrix = es_max,
    cluster_column = cluster_col,
    min_score_ratio = 0.25,
    return_all_scores = TRUE,
    top_n_per_cluster = TOP_N_CANDIDATES,
    verbose = FALSE
  )
  
  top_candidates <- sctype_results$all_scores %>%
    group_by(cluster) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    as.data.frame()
  
  cat(sprintf("  ✓ %d clusters, %d candidates\n\n", 
              length(unique(top_candidates$cluster)), nrow(top_candidates)))
  
  # 5.2 FindAllMarkers
  cat("Step 5.2: Finding cluster markers...\n")
  Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_col]]
  
  cluster_markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = LOGFC_THRESHOLD,
    test.use = "wilcox",
    verbose = FALSE
  )
  
  write.csv(cluster_markers, 
            file.path(res_dir, "cluster_markers_all.csv"), 
            row.names = FALSE)
  
  cluster_markers_sig <- cluster_markers %>%
    filter(p_val_adj < 0.05)
  
  write.csv(cluster_markers_sig, 
            file.path(res_dir, "cluster_markers_significant.csv"), 
            row.names = FALSE)
  
  cat(sprintf("  ✓ All markers: %d\n", nrow(cluster_markers)))
  cat(sprintf("  ✓ Significant markers: %d\n\n", nrow(cluster_markers_sig)))
  
  # 5.3 Compute marker overlaps for candidates
  cat("Step 5.3: Computing marker overlaps...\n")
  
  cluster_sig_genes <- cluster_markers_sig %>%
    group_by(cluster) %>%
    summarise(sig_genes = list(unique(gene)), .groups = "drop")
  
  candidates_with_overlap <- top_candidates %>%
    left_join(cluster_sig_genes, by = "cluster") %>%
    rowwise() %>%
    mutate(
      db_pos_genes = list(gs_list$gs_positive[[type]]),
      db_neg_genes = list(gs_list$gs_negative[[type]]),
      overlap_genes_pos = list(intersect(unlist(sig_genes), unlist(db_pos_genes))),
      overlap_genes_neg = list(intersect(unlist(sig_genes), unlist(db_neg_genes))),
      n_overlap_pos = length(overlap_genes_pos),
      n_overlap_neg = length(overlap_genes_neg),
      overlap_genes_pos_str = paste(overlap_genes_pos, collapse = ", "),
      overlap_genes_neg_str = paste(overlap_genes_neg, collapse = ", ")
    ) %>%
    ungroup() %>%
    as.data.frame()
  
  cat("  ✓ Marker overlaps computed\n\n")
  
  # 5.4 LLM validation
  cat("Step 5.4: LLM validation...\n")
  
  clusters <- unique(candidates_with_overlap$cluster)
  cat(sprintf("  Processing %d clusters sequentially\n", length(clusters)))
  
  serial_results <- lapply(clusters, function(cl) {
    result <- tryCatch({
      cand_data <- candidates_with_overlap %>%
        filter(cluster == cl) %>%
        mutate(
          celltype = type,
          overlap_genes_pos = overlap_genes_pos_str,
          overlap_genes_neg = overlap_genes_neg_str
        ) %>%
        select(rank, celltype, overlap_genes_pos, overlap_genes_neg)
      
      n_sig_markers <- sum(cluster_markers_sig$cluster == cl)
      
      prompt <- build_annotation_prompt(
        cluster_id = cl,
        n_cluster_markers = n_sig_markers,
        candidates_data = cand_data
      )
      
      response_text <- call_gpt_api(prompt)
      
      candidate_celltypes <- cand_data$celltype
      parsed_response <- parse_llm_response(response_text, candidate_celltypes)
      
      list(
        success = TRUE,
        cluster = cl,
        llm_result = data.frame(
          cluster = cl,
          selected_celltype = parsed_response$selected_celltype,
          confidence = parsed_response$confidence,
          reasoning = parsed_response$reasoning,
          key_markers_validated = paste(parsed_response$key_markers_validated, collapse = ", "),
          stringsAsFactors = FALSE
        ),
        llm_prompt = data.frame(
          cluster = cl,
          prompt = prompt,
          response = response_text,
          stringsAsFactors = FALSE
        )
      )
    }, error = function(e) {
      list(
        success = FALSE,
        cluster = cl,
        error_message = e$message
      )
    })
    return(result)
  })
  
  # Check for errors
  failed_clusters <- sapply(serial_results, function(x) !x$success)
  if (any(failed_clusters)) {
    cat("\n  ✗ ERROR: LLM API failed for the following clusters:\n")
    for (i in which(failed_clusters)) {
      cat(sprintf("    - Cluster %s: %s\n", 
                  serial_results[[i]]$cluster, 
                  serial_results[[i]]$error_message))
    }
    quit(save = "no", status = 1)
  }
  
  # Extract results
  llm_results_list <- lapply(serial_results, function(x) x$llm_result)
  llm_prompts_list <- lapply(serial_results, function(x) x$llm_prompt)
  
  # Combine results
  llm_decisions <- do.call(rbind, llm_results_list)
  llm_prompts <- do.call(rbind, llm_prompts_list)
  
  # Sort by cluster to ensure consistent ordering
  llm_decisions <- llm_decisions[order(llm_decisions$cluster), ]
  llm_prompts <- llm_prompts[order(llm_prompts$cluster), ]
  
  # Print completion summary
  cat(sprintf("\n  ✓ Successfully annotated %d clusters:\n", nrow(llm_decisions)))
  for (i in seq_len(nrow(llm_decisions))) {
    cat(sprintf("    - Cluster %s: %s\n", 
                llm_decisions$cluster[i], 
                llm_decisions$selected_celltype[i]))
  }
  
  # Attach overlap genes (pos/neg/all) for the selected celltype per cluster
  selected_overlaps <- candidates_with_overlap %>%
    mutate(celltype = type) %>%
    select(cluster, celltype, overlap_genes_pos_str, overlap_genes_neg_str) %>%
    distinct(cluster, celltype, .keep_all = TRUE)
  
  llm_decisions <- llm_decisions %>%
    left_join(selected_overlaps, by = c("cluster", "selected_celltype" = "celltype")) %>%
    mutate(
      overlap_genes_pos = ifelse(is.na(overlap_genes_pos_str), "", overlap_genes_pos_str),
      overlap_genes_neg = ifelse(is.na(overlap_genes_neg_str), "", overlap_genes_neg_str)
    ) %>%
    select(-overlap_genes_pos_str, -overlap_genes_neg_str) %>%
    rowwise() %>%
    mutate(
      overlap_genes_all = {
        pos_vec <- if (nzchar(overlap_genes_pos)) strsplit(overlap_genes_pos, ",\\s*")[[1]] else character(0)
        neg_vec <- if (nzchar(overlap_genes_neg)) strsplit(overlap_genes_neg, ",\\s*")[[1]] else character(0)
        all_vec <- unique(c(pos_vec, neg_vec))
        paste(all_vec, collapse = ", ")
      }
    ) %>%
    ungroup() %>%
    select(cluster, selected_celltype, confidence, reasoning, key_markers_validated,
           overlap_genes_all, overlap_genes_pos, overlap_genes_neg)
  
  write.csv(llm_decisions, 
            file.path(res_dir, "llm_decisions.csv"), 
            row.names = FALSE)
  write.csv(llm_prompts, 
            file.path(res_dir, "llm_prompts_log.csv"), 
            row.names = FALSE)
  
  cat("\n  ✓ LLM validation completed\n\n")
  
  # 5.5 Add annotations to Seurat object
  cat("Step 5.5: Adding annotations to Seurat object...\n")
  
  annotation_col <- paste0("llm_annotation_res", res)
  seurat_obj@meta.data[[annotation_col]] <- ""
  
  for (i in seq_len(nrow(llm_decisions))) {
    cl <- llm_decisions$cluster[i]
    celltype <- llm_decisions$selected_celltype[i]
    idx <- which(seurat_obj@meta.data[[cluster_col]] == cl)
    seurat_obj@meta.data[[annotation_col]][idx] <- celltype
  }
  
  cat("  ✓ Annotations added\n\n")
  
  # 5.6 Generate visualizations
  cat("Step 5.6: Generating visualizations...\n")
  
  # DimPlots
  plot_cluster_and_annotation(
    seurat_obj = seurat_obj,
    cluster_col = cluster_col,
    annotation_col = annotation_col,
    output_prefix = file.path(res_dir, "umap"),
    width = 14,
    height = 6
  )
  
  # DotPlot
  plot_marker_dotplot(
    seurat_obj = seurat_obj,
    llm_results = llm_decisions,
    cluster_col = cluster_col,
    output_prefix = file.path(res_dir, "umap"),
    max_markers_per_celltype = 10
  )
  
  cat("\n  ✓ Visualizations completed\n\n")
  
  cat(sprintf("  ✓ Resolution %s completed\n\n", res_str))
}

# ============================================================================
# Step 6: Save Final Results
# ============================================================================

cat("========================================================================\n")
cat("  Step 6: Saving Final Results\n")
cat("========================================================================\n\n")

cat("No Saving annotated Seurat object...\n")
# saveRDS(seurat_obj, file.path(opt$output_dir, "annotated_seurat.RDS"))
cat("  ✓ No Saved: annotated_seurat.RDS\n\n")

# Generate summary
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "mins")

summary_text <- paste0(
  "========================================================================\n",
  "  scMarkerAgent - Run Summary\n",
  "========================================================================\n\n",
  "Run Information:\n",
  sprintf("  Start time: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")),
  sprintf("  End time: %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")),
  sprintf("  Elapsed time: %.2f minutes\n\n", as.numeric(elapsed_time)),
  "Input:\n",
  sprintf("  File: %s\n", opt$input),
  sprintf("  Cells: %d\n", ncol(seurat_obj)),
  sprintf("  Features: %d\n\n", nrow(seurat_obj)),
  "Configuration:\n",
  sprintf("  Species: %s\n", opt$species),
  sprintf("  Tissues: %s\n", paste(tissues, collapse = ", ")),
  sprintf("  Disease types: %s\n", paste(disease_types, collapse = ", ")),
  sprintf("  Resolutions: %s\n", paste(resolution_strs, collapse = ", ")),
  sprintf("  Variable features: %d\n", opt$n_variable_features),
  sprintf("  PCA dimensions: 1:%d\n\n", opt$dims),
  "Database:\n",
  sprintf("  Cell types in gene sets: %d\n\n", length(gs_list$gs_positive)),
  "Results:\n",
  sprintf("  Output directory: %s\n", opt$output_dir),
  sprintf("  Resolutions processed: %d\n\n", length(resolutions)),
  "========================================================================\n"
)

writeLines(summary_text, file.path(opt$output_dir, "summary.txt"))
cat(summary_text)

cat("\n")
cat("========================================================================\n")
cat("  Pipeline Completed Successfully!\n")
cat("========================================================================\n\n")
cat(sprintf("Results saved to: %s\n\n", opt$output_dir))


