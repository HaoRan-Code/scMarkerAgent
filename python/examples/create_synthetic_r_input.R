#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})
args <- commandArgs(trailingOnly=FALSE)
match <- grep("^--file=", args, value=TRUE)
root <- dirname(normalizePath(sub("^--file=", "", match[1])))
counts <- readMM(file.path(root, "synthetic_counts.mtx"))
genes <- readLines(file.path(root, "synthetic_genes.tsv"), encoding="UTF-8")
cells <- readLines(file.path(root, "synthetic_cells.tsv"), encoding="UTF-8")
clusters <- read.delim(
  file.path(root, "synthetic_clusters.tsv"),
  stringsAsFactors=FALSE,
  check.names=FALSE
)
rownames(counts) <- genes
colnames(counts) <- cells
object <- CreateSeuratObject(counts=counts, project="scmarkeragent_synthetic")
object$cluster <- as.character(clusters$cluster[
  match(colnames(object), clusters$cell)
])
output <- file.path(root, "synthetic_input.rds")
saveRDS(object, output)
cat(output, "\n")
