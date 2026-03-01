> **Note**
> Large data files are not included in this GitHub repository.
> Please download the dataset from Zenodo: https://doi.org/10.5281/zenodo.18825921

This directory contains a **minimal set of input files** used by the example workflows. All files are RDS objects that can be loaded with R `readRDS()`.

## File descriptions

- **`scmarkeragent_db.RDS`**
  
  - **What it is**: The scMarkerAgent marker annotation resource used by the scripts to retrieve positive/negative marker gene sets after filtering by `species`, `tissue`, `disease_name`, `cell_name`, and `polarity`.
  - **Key content**: Built by mining **294,692** PubMed Central full-text articles with a semi-automated LLM-driven framework, extracting **1,352,754** candidate cell type–marker statements with sentence-level evidence and applying evidence-based quality control; after filtering to human/mouse/rat, **890,296** high-quality annotations were retained, including **82,165** curated negative marker annotations and **417,812** disease-context annotations. Tissues/cell types/diseases are standardized with Uberon/Cell Ontology/Disease Ontology, and marker names are standardized via NCBI E-utilities (with additional species-specific resources and UniProt reviewed entries).
  - **Used by**:
    - `cell_annotation/annotate_cells.R`: prepares marker gene sets for ScType scoring and subsequent LLM-based validation
    - `cell_score/ucell_demo.R`: builds a marker signature from the filtered database and computes per-cell scores with UCell

- **`example_pbmc_input.rds`**
  
  - **What it is**: An example Seurat object (PBMC) used as input for the **Cell Annotation** module.
  - **Used by**: `cell_annotation/annotate_cells.R --input`
  - **Reference**: Zheng, Grace X Y et al. “Massively parallel digital transcriptional profiling of single cells.” *Nature Communications* 8 (2017): 14049. doi:10.1038/ncomms14049.

- **`example_glioblastoma_input.rds`**
  
  - **What it is**: An example Seurat object (glioblastoma) used as input for the **Cell Scoring** module.
  - **Used by**: `cell_score/score_cells.R --input_seurat`
  - **Reference**: Ruiz-Moreno, Cristian et al. “Charting the single-cell and spatial landscape of IDH-wild-type glioblastoma with GBmap.” *Neuro-Oncology* 27.9 (2025): 2281–2295. doi:10.1093/neuonc/noaf113.
