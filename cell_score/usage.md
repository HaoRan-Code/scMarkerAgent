# Cell Scoring Module

Per-cell marker gene scoring using the UCell algorithm, powered by the scMarkerAgent database.

## Usage Examples

```bash
# Negative marker scoring — 4 threads
Rscript ./cell_score/ucell_demo.R \
  --species Human \
  --tissue brain \
  --disease_name glioblastoma \
  --cell_name "malignant cell" \
  --polarity negative \
  --n_threads 4 \
  --outdir ./results/cell_score_negative_example_4threads \
  --input_seurat ./scMarkerAgent_data/example_glioblastoma_input.rds \
  --marker_rds ./scMarkerAgent_data/scmarkeragent_db.RDS

# Negative marker scoring — 6 threads
Rscript ./cell_score/ucell_demo.R \
  --species Human \
  --tissue brain \
  --disease_name glioblastoma \
  --cell_name "malignant cell" \
  --polarity negative \
  --n_threads 6 \
  --outdir ./results/cell_score_negative_example_6threads \
  --input_seurat ./scMarkerAgent_data/example_glioblastoma_input.rds \
  --marker_rds ./scMarkerAgent_data/scmarkeragent_db.RDS

# Positive marker scoring — 4 threads
Rscript ./cell_score/ucell_demo.R \
  --species Human \
  --tissue brain \
  --disease_name glioblastoma \
  --cell_name "malignant cell" \
  --polarity positive \
  --n_threads 4 \
  --outdir ./results/cell_score_positive_example_4threads \
  --input_seurat ./scMarkerAgent_data/example_glioblastoma_input.rds \
  --marker_rds ./scMarkerAgent_data/scmarkeragent_db.RDS

# Positive marker scoring — 6 threads
Rscript ./cell_score/ucell_demo.R \
  --species Human \
  --tissue brain \
  --disease_name glioblastoma \
  --cell_name "malignant cell" \
  --polarity positive \
  --n_threads 6 \
  --outdir ./results/cell_score_positive_example_6threads \
  --input_seurat ./scMarkerAgent_data/example_glioblastoma_input.rds \
  --marker_rds ./scMarkerAgent_data/scmarkeragent_db.RDS
```
