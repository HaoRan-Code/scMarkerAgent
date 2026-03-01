# Cell Annotation Agent

Automated cell type annotation pipeline integrating the scMarkerAgent database with LLM-based validation.

## Quick Start

```bash
# View help
Rscript ./cell_annotation/annotate_cells.R --help
Rscript ./cell_annotation/annotate_cells.R --help_extended

# Basic usage
Rscript ./cell_annotation/annotate_cells.R \
  --input ./scMarkerAgent_data/example_pbmc_input.rds \
  --tissue blood \
  --output_dir ./results/my_annotation
```

## Full Usage Examples

```bash
# Disease + Normal
Rscript ./cell_annotation/annotate_cells.R \
  --input ./scMarkerAgent_data/example_pbmc_input.rds \
  --species Human \
  --tissue "blood" \
  --n_variable_features 2000 \
  --dims 30 \
  --cluster_resolutions "0.3,0.5,1.0" \
  --min_markers_per_cell 2 \
  --disease_type "Normal,acute myeloid leukemia,aplastic anemia" \
  --output_dir ./results/example_full \
  --random_seed 1234

# Normal only
Rscript ./cell_annotation/annotate_cells.R \
  --input ./scMarkerAgent_data/example_pbmc_input.rds \
  --species Human \
  --tissue "blood" \
  --n_variable_features 2000 \
  --dims 30 \
  --cluster_resolutions "0.3,0.5,1.0" \
  --min_markers_per_cell 2 \
  --disease_type "Normal" \
  --output_dir ./results/example_normal \
  --random_seed 1234

# Disease only
Rscript ./cell_annotation/annotate_cells.R \
  --input ./scMarkerAgent_data/example_pbmc_input.rds \
  --species Human \
  --tissue "blood" \
  --n_variable_features 2000 \
  --dims 30 \
  --cluster_resolutions "0.3,0.5,1.0" \
  --min_markers_per_cell 2 \
  --disease_type "acute myeloid leukemia" \
  --output_dir ./results/example_disease \
  --random_seed 1234
```

## Parameters

### Required
| Parameter | Description |
|-----------|-------------|
| `--input` / `-i` | Path to input Seurat object (RDS file) |
| `--tissue` | Tissue type (e.g., blood, brain, lung) |

### scRNA-seq Processing
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--n_variable_features` | 2000 | Number of highly variable features |
| `--dims` | 30 | Number of PCA dimensions (1:N) |
| `--cluster_resolutions` | 0.5 | Clustering resolutions, comma-separated (max 3) |
| `--random_seed` | 1234 | Random seed for reproducibility |

### Database Filtering
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--species` | Human | Species: Human, Mouse, or Rat |
| `--min_markers_per_cell` | 2 | Minimum number of markers per cell type |
| `--disease_type` | Normal | Disease types, comma-separated |

### Output
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output_dir` / `-o` | ./results | Output directory |

### Fixed Internal Parameters
| Parameter | Value |
|-----------|-------|
| `logfc_threshold` | 1 |
| `top_n_candidates` | 10 |
| `assay` | RNA |
| `k_param` | 20 |

## Output Structure

```
results/
├── summary.txt                       # Run summary
└── resolution_X.X/                   # Results per resolution
    ├── cluster_markers_all.csv       # All cluster markers
    ├── cluster_markers_significant.csv  # Significant markers (adj. p < 0.05)
    ├── llm_decisions.csv             # LLM annotation decisions
    ├── llm_prompts_log.csv           # Full LLM interaction log
    ├── umap_cluster_annotation.pdf   # UMAP visualization
    └── umap_marker_dotplot.pdf       # Marker dot plot
```

## API Configuration

Before running, set your OpenAI API key in `annotate_cells.R`:

```r
Sys.setenv(PRIMARY_API_KEY = "YOUR_API_KEY_HERE")
Sys.setenv(PRIMARY_API_URL = "https://api.openai.com/v1/chat/completions")
Sys.setenv(PRIMARY_MODEL = "gpt-4o")
```

## Backend Integration

```python
import subprocess

def annotate_cells(input_path, tissue, output_dir, **kwargs):
    cmd = [
        "Rscript",
        "./cell_annotation/annotate_cells.R",
        "--input", input_path,
        "--tissue", tissue,
        "--output_dir", output_dir
    ]
    for key, value in kwargs.items():
        cmd.extend([f"--{key}", str(value)])
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise Exception(f"Annotation failed: {result.stderr}")
    return output_dir

# Example
annotate_cells(
    input_path="./scMarkerAgent_data/example_pbmc_input.rds",
    tissue="blood",
    output_dir="./results/pbmc"
)
```

## Notes

- Input file size limit: ≤ 5 GB
- Requires network access for LLM API calls
- LLM API failure will terminate the pipeline

## Module Structure

```
cell_annotation/
├── annotate_cells.R          # Main script
└── functions/
    ├── sctype_improved_functions.R   # ScType scoring algorithm
    ├── data_preparation.R            # Data preparation utilities
    ├── llm_api.R                     # LLM API interface
    └── visualization.R               # UMAP and dot plot generators
```
