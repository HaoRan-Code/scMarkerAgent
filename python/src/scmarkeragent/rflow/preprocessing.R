#!/usr/bin/env Rscript
# =============================================================================
# preprocessing.R -- preprocessing for the .rds arm.
#
# Both arms START FROM RAW COUNTS and run the SAME standard single-cell
# workflow (so the .h5ad and .rds pipelines are methodologically identical
# -- no reuse of any pre-computed normalization / clustering from the input
# object):
#   1. take RAW counts from the EXPLICITLY DECLARED assay/layer (SCMA_COUNTS_SOURCE)
#   2. sum case-insensitive duplicate gene symbols
#   3. QC: CreateSeuratObject(min.cells=3, min.features=200) + percent.mt < 20
#   4. normalize: LogNormalize, scale.factor 1e4   (== CP10k + log1p)
#   5. cluster: 2000 vst HVG -> ScaleData -> PCA(30) -> SNN kNN(k=20)
#              -> Leiden (FindClusters algorithm=4, leiden_method="igraph", res).
#      The partition is ALWAYS recomputed; a partition carried by the input object is
#      never reused.
#   6. optional UMAP controlled by the shared configuration
#   7. per-cluster DE: one genome-wide presto::wilcoxauc pass on the freshly
#      log-normalized data over ALL cells -- the gene universe and per-cluster BH family
#      every arm uses. The menu table is a slice of that pass over
#      DB-menu(positive + negative) INTERSECT measured genes. avg_log2FC is added with the
#      Seurat-style formula both arms use, because it is what the gate reads.
#
# The shared random state drives Leiden and dimensional reduction.
#
# Outputs (cache/, consumed by candidate_scoring.R and the shared report layer):
#   <tag>_de_meta.rds : list(de, meta(cell,cluster),
#                            genes, cluster_col, species, tissue, disease, umap)
#   <tag>_emb.rds     : PCA embedding (cells x <=20), rownames = cells
#
#   Usage: Rscript preprocessing.R <tag>
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})
.script_argument <- grep(
  "^--file=", commandArgs(FALSE),
  value = TRUE
)
if (length(.script_argument) != 1L) {
  stop("preprocessing.R must be executed with Rscript")
}
.sd <- dirname(normalizePath(sub("^--file=", "", .script_argument)))
source(file.path(.sd, "config.R"))
source(file.path(.sd, "marker_database.R"))

# ---- generic production input contract -------------------------------------
tag <- INPUT_TAG
INPUT_RDS <- .check_input(INPUT_RDS, "r")
if (!nzchar(INPUT_SPECIES) || !nzchar(INPUT_TISSUE)) {
  stop("SCMA_SPECIES and SCMA_TISSUE are required")
}
if (!nzchar(COUNTS_SOURCE)) {
  stop("SCMA_COUNTS_SOURCE is required, as '<assay>/<layer>' (e.g. 'RNA/counts')")
}
cfg <- list(
  obj = INPUT_RDS,
  species = INPUT_SPECIES, tissue = INPUT_TISSUE, disease = INPUT_DISEASE,
  res = LEIDEN_RESOLUTION
)
RANDOM_STATE <- as.integer(CFG$preprocessing$random_state)
COMPUTE_UMAP <- isTRUE(CFG$features$compute_umap)
PCA_N <- as.integer(CFG$preprocessing$pca_n)
NEIGHBORS_K <- as.integer(CFG$preprocessing$neighbors_k)
RES <- as.numeric(cfg$res)
set.seed(RANDOM_STATE)

cat(sprintf(
  "==== prep(rds): %s  (%s/%s/%s) counts_source=%s random_state=%d res=%g umap=%s ====\n",
  tag, cfg$species, cfg$tissue, paste(cfg$disease, collapse = "+"),
  COUNTS_SOURCE, RANDOM_STATE, RES, COMPUTE_UMAP
))

# ---- raw counts from the EXPLICITLY DECLARED assay/layer -------------------
# Nothing is inferred. The declaration is echoed into the run configuration, so a
# normalized matrix can never be silently annotated as if it were counts.
.counts_parts <- strsplit(COUNTS_SOURCE, "/", fixed = TRUE)[[1]]
if (length(.counts_parts) != 2L || !all(nzchar(.counts_parts))) {
  stop("SCMA_COUNTS_SOURCE must be '<assay>/<layer>', e.g. 'RNA/counts'; got: ", COUNTS_SOURCE)
}
.counts_assay <- .counts_parts[1]
.counts_layer <- .counts_parts[2]
o0 <- readRDS(cfg$obj)
if (!(.counts_assay %in% names(o0@assays))) {
  stop(sprintf(
    "counts source '%s' requested but assay '%s' is absent; available assays: %s",
    COUNTS_SOURCE, .counts_assay, paste(names(o0@assays), collapse = ", ")
  ))
}
DefaultAssay(o0) <- .counts_assay
cat("loaded", paste(dim(o0), collapse = " x "), "(genes x cells)\n")
counts <- tryCatch(
  GetAssayData(o0, assay = .counts_assay, layer = .counts_layer),
  error = function(e) NULL
)
if (is.null(counts) || !nrow(counts)) {
  stop("counts source '", COUNTS_SOURCE, "' is absent or empty in object: ", cfg$obj)
}
.count_values <- if (inherits(counts, "sparseMatrix")) counts@x else as.numeric(counts)
if (any(!is.finite(.count_values)) || any(.count_values < 0)) {
  stop("counts source '", COUNTS_SOURCE, "' contains non-finite or negative values: ", cfg$obj)
}
if (length(.count_values) && max(abs(.count_values - round(.count_values))) > 1e-6) {
  stop(
    "counts source '", COUNTS_SOURCE, "' is not integer-like; point SCMA_COUNTS_SOURCE at ",
    "the raw counts, or round the input once at its origin: ", cfg$obj
  )
}
rm(.count_values)
md0 <- o0@meta.data
rm(o0)
gc()

# ---- sum CASE-INSENSITIVE duplicate gene symbols ---------------------------
# A source that carries the same gene under two casings (e.g. both Pisd and PISD) must
# contribute ONE summed row, or the two arms disagree on that gene's counts and on the
# measured-gene universe. The NATIVE symbol (first seen) is kept for storage, so species
# symbol nomenclature is never mutated.
.sym <- rownames(counts)
.key <- toupper(.sym)
if (anyDuplicated(.key)) {
  .ukey <- unique(.key)
  .collapse <- Matrix::sparseMatrix(
    i = match(.key, .ukey), j = seq_along(.key),
    x = rep(1, length(.key)), dims = c(length(.ukey), length(.key))
  )
  counts <- .collapse %*% counts
  rownames(counts) <- .sym[match(.ukey, .key)] # native representative per upper key
  cat(sprintf(
    "  gene symbols: summed %d case-insensitive duplicates -> %d unique genes\n",
    length(.key) - length(.ukey), length(.ukey)
  ))
  rm(.ukey, .collapse)
}
rm(.sym, .key)

# ---- standard Seurat workflow FROM RAW -------------------------------------
n0 <- ncol(counts)
o <- CreateSeuratObject(counts = counts, meta.data = md0, min.cells = QC_MIN_CELLS, min.features = QC_MIN_GENES)
mito_feats <- rownames(o)[grepl("^MT-", toupper(rownames(o)))] # MT- prefix on uppercased symbols
o[["percent.mt"]] <- if (length(mito_feats)) PercentageFeatureSet(o, features = mito_feats) else 0
o <- subset(o, subset = percent.mt < QC_MAX_MT)
cat(sprintf(
  "  QC: %d -> %d cells, %d genes (min_genes=%d, min_cells=%d, max_mt=%g)\n",
  n0, ncol(o), nrow(o), QC_MIN_GENES, QC_MIN_CELLS, QC_MAX_MT
))
o <- NormalizeData(o,
  normalization.method = "LogNormalize",
  scale.factor = as.numeric(CFG$preprocessing$normalization_scale_factor),
  verbose = FALSE
)
o <- FindVariableFeatures(
  o,
  selection.method = "vst", nfeatures = as.integer(CFG$preprocessing$hvg_n), verbose = FALSE
)
o <- ScaleData(o, verbose = FALSE)
o <- RunPCA(
  o,
  npcs = PCA_N, verbose = FALSE, seed.use = RANDOM_STATE
)
# Clustering is ALWAYS recomputed; a partition carried by the input object is never
# reused, so every dataset and both arms enter scoring from the same workflow.
o <- FindNeighbors(o, dims = seq_len(PCA_N), k.param = NEIGHBORS_K, verbose = FALSE)
set.seed(RANDOM_STATE + 1L)
o <- FindClusters(
  o,
  resolution = RES, algorithm = 4, leiden_method = "igraph",
  random.seed = RANDOM_STATE + 1L, verbose = FALSE
)
clusters <- as.character(o$seurat_clusters)
cat(sprintf(
  "  clusters: Leiden(igraph) res=%g -> %d clusters\n",
  RES, length(unique(clusters))
))
cluster_col_out <- sprintf("leiden_res%g", RES)
umap <- NULL
if (COMPUTE_UMAP) {
  o <- RunUMAP(
    o,
    dims = seq_len(PCA_N), seed.use = RANDOM_STATE, verbose = FALSE
  )
  umap <- Embeddings(o, "umap")
}
emb <- Embeddings(o, "pca")[, seq_len(min(
  as.integer(CFG$preprocessing$score_pca_n),
  ncol(Embeddings(o, "pca"))
)), drop = FALSE]

# ---- DE gene universe = DB-menu (positive markers) INTERSECT measured genes -----------------
# The shared gene universe, narrowed by the candidate_scoring.R database context filter.
# Both arms compute per-cluster DE over the SAME menu-measured family, so the per-cluster
# BH denominator matches. Testing all measured genes here instead would give this arm a
# far larger and systematically more conservative BH family, breaking parity.
# menu_pos = positive-marker gene symbols for this (species, tissue, disease) context (NOT the
# eligibility-gated scoring shortlist -- that gating happens later in candidate_scoring.R); the DE
# universe is therefore a superset of every gene the scorer looks up. cross_species pooling
# Cross-species mode from the shared configuration is honored so the DE universe still
# covers the pooled ortholog genes; off by default.
.menu_measured_genes <- function(measured, sp_q, ti_q, di_q) {
  suppressPackageStartupMessages({
    library(data.table)
  })
  source(file.path(SELF, "uberon_ontology.R"))
  X_SPECIES <- isTRUE(CFG$features$cross_species_markers)
  INVITRO_PAT <- paste0(
    "induced pluripotent stem cell|pluripotent stem cell|embryonic stem cell|",
    "\\bipsc\\b|\\bipscs\\b|\\bhipsc\\b|\\bips cell\\b|\\besc\\b|\\bcell line\\b|\\bcell-line\\b|",
    "\\borganoid\\b|\\bspheroid\\b|\\bteratoma\\b|\\bin vitro\\b|\\bin-vitro\\b|t-cell line"
  )
  XSP <- c(
    Human = "\\bmouse\\b|\\bmurine\\b|mus musculus|\\brat\\b|\\brats\\b|rattus|zebrafish|drosophila|xenopus|\\bporcine\\b|\\bbovine\\b|\\bcanine\\b|\\bovine\\b",
    Mouse = "\\bhuman\\b|\\brat\\b|\\brats\\b|rattus|zebrafish|drosophila|xenopus",
    Rat = "\\bhuman\\b|\\bmouse\\b|\\bmurine\\b|mus musculus|zebrafish|drosophila|xenopus"
  )
  HOMOLOGY_PAT <- "equivalent of|homolog|ortholog|counterpart"
  UB <- uberon_load(OBO_UBERON)
  tissue_closure <- function(ti, tissue_root = TISSUE_ROOT) { # == candidate_scoring.R .tissue_closure
    tr <- tolower(trimws(if (!nzchar(tissue_root)) "self" else tissue_root))
    ids <- character(0)
    nms <- .ub_norm(ti)
    rid <- uberon_resolve(UB, ti)
    if (!is.na(rid)) ids <- rid
    if (tr %in% c("exact", "none", "off")) {
      return(list(ids = ids, names = nms))
    }
    root <- if (tr %in% c("self", "")) ti else tissue_root
    cids <- uberon_closure_ids(UB, root)
    if (length(cids)) {
      ids <- unique(c(ids, cids))
      nms <- unique(c(nms, uberon_closure_names(UB, root)))
    }
    list(ids = ids, names = nms)
  }
  db <- scma_load_marker_db(DB_CSV)
  setDT(db)
  .ctl <- tolower(db$cell_type)
  .mm <- rep(FALSE, nrow(db)) # task4 cross-species-mismatch source cleanup
  for (.sp in names(XSP)) {
    .i <- which(db$species == .sp)
    .mm[.i] <- grepl(XSP[[.sp]], .ctl[.i], perl = TRUE) & !grepl(HOMOLOGY_PAT, .ctl[.i], perl = TRUE)
  }
  if (any(.mm)) db <- db[!.mm]
  if ("is_in_vitro" %in% names(db)) db[, is_in_vitro := tolower(as.character(is_in_vitro)) %in% c("true", "1", "yes", "t")] else db[, is_in_vitro := grepl(INVITRO_PAT, tolower(cell_type), perl = TRUE)]
  db[, .ttn := .ub_norm(tissue_type)]
  .TC <- tissue_closure(ti_q)
  tt_match <- (db$tissue_type_uberon_id %in% .TC$ids) | (db[[".ttn"]] %in% .TC$names)
  iv_ok <- if (EXCLUDE_IN_VITRO) !db$is_in_vitro else rep(TRUE, nrow(db))
  .dmatch <- .disease_exact_match(db$disease_normalized, di_q)
  base_ok <- tt_match & iv_ok & .dmatch &
    (toupper(db$gene_qc_pass) %in% c("", "TRUE") | is.na(db$gene_qc_pass)) &
    !is.na(db$gene_symbol) & !(tolower(db$gene_symbol) %in% c("", "unknown")) &
    !is.na(db$cell_type) & db$cell_type != "" &
    db$marker_polarity %in% c("positive", "negative")
  menu_pos <- unique(as.character(db[species == sp_q & base_ok & marker_polarity == "positive", gene_symbol]))
  menu_neg <- unique(as.character(db[species == sp_q & base_ok & marker_polarity == "negative", gene_symbol]))
  if (X_SPECIES) { # pool OTHER species' orthologs (== candidate_scoring.R)
    source(file.path(SELF, "ortho_map.R"))
    for (osp in setdiff(c("Human", "Mouse", "Rat"), sp_q)) {
      pr <- ortho_pairs(osp, sp_q, ORTHO_DIR)
      if (!nrow(pr)) next
      og_pos <- unique(toupper(as.character(db[species == osp & base_ok & marker_polarity == "positive", gene_symbol])))
      og_neg <- unique(toupper(as.character(db[species == osp & base_ok & marker_polarity == "negative", gene_symbol])))
      if (length(og_pos)) menu_pos <- unique(c(menu_pos, pr[src_sym %in% og_pos, tgt_sym]))
      if (length(og_neg)) menu_neg <- unique(c(menu_neg, pr[src_sym %in% og_neg, tgt_sym]))
    }
  }
  meas <- as.character(measured)
  # CASE-SENSITIVE menu-vs-measured intersection. A gene symbol IS its
  # case: HGNC PECAM1 and MGI Pecam1 name two species' measurements, and matching them
  # on an uppercased key silently annotated a mouse matrix against the human menu. An
  # exact match is the only one that means the same gene, so the matched symbol is both
  # the curated symbol and the matrix's own native symbol.
  measured_native <- function(symbols) {
    sort(intersect(unique(as.character(symbols)), meas))
  }
  positive_genes <- measured_native(menu_pos)
  negative_genes <- setdiff(measured_native(menu_neg), positive_genes)
  if (!length(positive_genes)) {
    # The symbols from both sides, so the reader can see WHY nothing matched instead of
    # being told to re-check three dropdowns. Nothing here infers what species the object
    # is: it reports what was measured, what is curated, and what to do.
    .sample <- function(symbols, cap = 4L) {
      shown <- head(sort(unique(as.character(symbols))), cap)
      if (!length(shown)) "nothing" else paste(shown, collapse = ", ")
    }
    stop(sprintf(
      paste0(
        "not one of the %d curated positive markers for %s/%s/%s is among this object's ",
        "%d measured genes. Gene symbols are matched exactly, case included -- this ",
        "object measures %s; the menu curates %s. If those are the same genes written ",
        "differently, this object is not written in the nomenclature the selected ",
        "species is curated under: select the species these symbols belong to, or ",
        "rewrite the symbols in the selected species' nomenclature. Otherwise the ",
        "tissue or disease does not describe this object."
      ),
      length(unique(as.character(menu_pos))), sp_q, ti_q, paste(di_q, collapse = "+"),
      length(meas), .sample(meas), .sample(menu_pos)
    ), call. = FALSE)
  }
  cat(sprintf(
    "  menu genes: positive=%d negative=%d, measured=%d, menu-measured=%d (%d positive + %d negative-only)%s\n",
    length(unique(as.character(menu_pos))), length(unique(as.character(menu_neg))), length(meas),
    length(positive_genes) + length(negative_genes), length(positive_genes),
    length(negative_genes), if (X_SPECIES) " [cross_species ON]" else ""
  ))
  list(
    positive = positive_genes, negative = negative_genes,
    all = sort(unique(c(positive_genes, negative_genes)))
  )
}
# ---- per-cluster DE on the freshly normalized data (presto, DB-menu-measured genes) ------
# EVERY cell is tested. Subsampling made the DE table -- and therefore the significance
# gate -- depend on a random draw, which is not something an annotation should depend on.
# The freshly built object always carries one "RNA" assay (CreateSeuratObject default);
# .counts_assay named where the counts were READ FROM in the input object, not here.
mat <- GetAssayData(o, assay = "RNA", layer = "data") # full log-norm (genes x cells); KEPT for the norm cache below
# The scaled and raw-count layers exist only for HVG selection and PCA, which are done.
# Dropping them here is what keeps a 500k-cell input from holding three copies of the
# matrix through the two DE passes; `cell_names` is captured first because the object is
# no longer consulted for it afterwards.
cell_names <- colnames(o)
try(
  {
    o[["RNA"]]$scale.data <- NULL
    o[["RNA"]]$counts <- NULL
  },
  silent = TRUE
)
gc(verbose = FALSE)
.menu <- .menu_measured_genes(rownames(mat), cfg$species, cfg$tissue, cfg$disease)
de_genes <- .menu$all

# ---- ONE genome-wide one-vs-rest DE pass -----------------------------------------------
# Every gene is tested in a single pass, so the per-cluster BH family is the whole
# transcriptome and the menu table below is a slice of that same pass rather than a second
# presto run carrying its own, smaller family. presto::wilcoxauc returns `pval` (two-sided
# Mann-Whitney normal approximation with a tie-corrected rank variance) and a `padj` that is
# already per-group BH; it is recomputed explicitly so the semantics hold regardless of the
# presto version, and it is the column the retrieval gate reads.
cat(sprintf("  genome-wide DE (presto::wilcoxauc) on %d genes x %d cells ...\n", nrow(mat), ncol(mat)))
t1 <- Sys.time()
full_de <- presto::wilcoxauc(mat, clusters)
setDT(full_de)
full_de[, padj := p.adjust(pval, method = "BH"), by = group]
full_de[, group := as.character(group)]
full_de[, feature := as.character(feature)]
full_de[, m_min := min(avgExpr), by = feature]
full_de[, m_max := max(avgExpr), by = feature]
full_de[, scaled_mean := fifelse(m_max > m_min, (avgExpr - m_min) / (m_max - m_min), 0)]
full_de[, c("m_min", "m_max") := NULL]

# avg_log2FC on top of wilcoxauc: de-log1p first, compare arithmetic means,
# add a pseudocount of one, take log2. presto's own `logFC` is a difference of mean
# LOG-normalized expression and is kept for audit only; avg_log2FC is what the shared
# retrieval gate reads, so it has to be produced identically by both arms.
.linear <- mat
.linear@x <- expm1(.linear@x)
.total <- Matrix::rowSums(.linear)
.n_all <- ncol(.linear)
.lfc <- rbindlist(lapply(unique(clusters), function(cl) {
  index <- which(clusters == cl)
  sum_in <- Matrix::rowSums(.linear[, index, drop = FALSE])
  mean_in <- sum_in / length(index)
  n_out <- .n_all - length(index)
  mean_out <- if (n_out > 0L) (.total - sum_in) / n_out else mean_in
  data.table(
    feature = rownames(.linear), group = as.character(cl),
    avg_log2FC = as.numeric(log2((mean_in + 1.0) / (mean_out + 1.0)))
  )
}))
full_de <- merge(full_de, .lfc, by = c("feature", "group"), all.x = TRUE, sort = FALSE)
if (anyNA(full_de$avg_log2FC)) stop("avg_log2FC could not be resolved for every DE row")
rm(.linear, .total, .lfc)
gc(verbose = FALSE)
cat("  genome-wide DE rows:", nrow(full_de), "in",
  round(as.numeric(Sys.time() - t1, units = "secs"), 1), "s\n")

# menu slice, unfiltered, feature-major over de_genes (== preprocessing.menu_de_table)
de <- full_de[feature %in% de_genes]
de[, .feature_order := match(feature, de_genes)]
setorder(de, .feature_order, group)
de[, .feature_order := NULL]

# ---- FindAllMarkers-equivalent reporting screen ----------------------------------------
# Seurat's `return.thresh` is deliberately NOT applied: this is the "all markers" half of a
# pair whose companion is the gate-passing subset.
.screen <- CFG$preprocessing$all_markers
markers_all <- full_de[
  pmax(pct_in, pct_out) >= 100 * as.numeric(.screen$min_pct) &
    abs(avg_log2FC) >= as.numeric(.screen$logfc_threshold)
]
cat(sprintf(
  "  all-gene DE: %d reported rows (min_pct=%g, logfc_threshold=%g)\n",
  nrow(markers_all), as.numeric(.screen$min_pct), as.numeric(.screen$logfc_threshold)
))
rm(full_de)
gc(verbose = FALSE)

meta <- data.table(
  cell = cell_names, cluster = clusters
)

# NATIVE log-norm matrix (genes x cells) for candidate_scoring.R (so it need NOT reload
# the big object). Gene casing is PRESERVED (species
# symbol nomenclature); candidate_scoring.R matches DB genes by EXACT symbol.
dat_u <- mat[de_genes, , drop = FALSE]
saveRDS(
  list(
    de = de, meta = meta, genes = rownames(mat), cluster_col = cluster_col_out,
    menu_genes = de_genes,
    menu_positive_genes = .menu$positive, menu_negative_genes = .menu$negative,
    leiden_resolution = RES, counts_source = COUNTS_SOURCE,
    source_path = normalizePath(cfg$obj),
    species = cfg$species, tissue = cfg$tissue, disease = cfg$disease, umap = umap
  ),
  file.path(OUT, paste0(tag, "_de_meta.rds"))
)
saveRDS(emb, file.path(OUT, paste0(tag, "_emb.rds")))
saveRDS(dat_u, file.path(OUT, paste0(tag, "_norm.rds")))
saveRDS(markers_all, file.path(OUT, paste0(tag, "_markers_all.rds")))
cat(
  "[done] ->", paste0(tag, "_de_meta.rds"), "+", paste0(tag, "_emb.rds"), "+",
  paste0(tag, "_markers_all.rds"), "| clusters:", length(unique(meta$cluster)),
  "| all-gene rows:", nrow(markers_all), "\n"
)
