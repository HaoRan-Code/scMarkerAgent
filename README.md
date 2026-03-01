> **Note**
> Large data files are not included in this GitHub repository.
> Data can be downloaded from Zenodo: https://doi.org/10.5281/zenodo.18825921

# scMarkerAgent: An LLM Evidence Agent-based Cell Marker Atlas

**scMarkerAgent** is an evidence-based cell marker intelligence agent built through an LLM-driven literature curation framework for transparent and high-specificity single-cell RNA-seq cell type annotation.

## Overview

The scMarkerAgent database was constructed by mining **294,692** PubMed Central full-text articles using a semi-automated multi-agent framework, extracting **1,352,754** candidate cell type–marker statements with sentence-level evidence. After evidence-based quality control, **890,296** high-quality annotations were retained, including **82,165** curated negative marker annotations and **417,812** disease-context annotations. Tissues, cell types, and diseases are standardized against Uberon, Cell Ontology, and Disease Ontology; marker gene symbols are normalized via NCBI E-utilities, species-specific resources, and UniProt reviewed entries.

Every annotation is traceable to sentence-level literature evidence and can be audited in biological context. The resource is released as a FAIR-compliant database with a code-free web platform for marker retrieval, automated cell annotation, and customizable cell scoring.

## Repository Structure

```
scMarkerAgent/
├── README.md                          # This file
├── cell_annotation/                   # Cell Annotation Agent module
│   ├── usage.md
│   ├── annotate_cells.R               # Main annotation pipeline script
│   └── functions/
│       ├── sctype_improved_functions.R
│       ├── data_preparation.R
│       ├── llm_api.R
│       └── visualization.R
├── cell_score/                        # Cell Scoring module
│   ├── usage.md
│   └── ucell_demo.R                   # UCell-based per-cell scoring script
└── scMarkerAgent_data/                # Input data files
    └── README.md
```

## Modules

### Cell Annotation Agent (`cell_annotation/`)

An end-to-end automated pipeline for cell type annotation of scRNA-seq data. The workflow integrates:

1. Standard Seurat preprocessing (normalization, PCA, clustering, UMAP)
2. scMarkerAgent database query to retrieve positive and negative marker gene sets
3. ScType scoring for candidate cell type ranking per cluster
4. LLM-based validation to resolve ambiguities and produce final annotations with confidence scores and supporting evidence

See [`cell_annotation/usage.md`](cell_annotation/usage.md) for full documentation.

**Quick start:**

```bash
Rscript ./cell_annotation/annotate_cells.R \
  --input ./scMarkerAgent_data/example_pbmc_input.rds \
  --tissue blood \
  --output_dir ./results/pbmc_annotation
```

### Cell Scoring (`cell_score/`)

A standalone UCell-based per-cell scoring script that retrieves a marker signature from the scMarkerAgent database and computes enrichment scores at single-cell resolution. Supports both positive and negative marker polarities.

See [`cell_score/usage.md`](cell_score/usage.md) for usage examples.

## Web Platform

The scMarkerAgent web platform is available at:

- http://www.markeragent.net

## Data

Input data files are stored in `scMarkerAgent_data/`. See [`scMarkerAgent_data/README.md`](scMarkerAgent_data/README.md) for descriptions of each file.

| File                             | Description                                                    |
| -------------------------------- | -------------------------------------------------------------- |
| `scmarkeragent_db.RDS`           | scMarkerAgent marker annotation database (890,296 annotations) |
| `example_pbmc_input.rds`         | Example PBMC Seurat object for cell annotation                 |
| `example_glioblastoma_input.rds` | Example glioblastoma Seurat object for cell scoring            |

## API Configuration

The Cell Annotation Agent requires an OpenAI-compatible API key. Before running, set the key in `cell_annotation/annotate_cells.R`:

```r
Sys.setenv(PRIMARY_API_KEY = "YOUR_API_KEY_HERE")
Sys.setenv(PRIMARY_API_URL = "https://api.openai.com/v1/chat/completions")
Sys.setenv(PRIMARY_MODEL = "gpt-4o")
```

## Dependencies

### R Packages

- **Seurat** — single-cell analysis
- **UCell** — per-cell scoring
- **optparse** — command-line argument parsing
- **httr**, **jsonlite** — LLM API communication
- **ggplot2**, **patchwork**, **ggrastr** — visualization
- **dplyr**, **data.table** — data manipulation
- **scales**, **showtext** — plot utilities
