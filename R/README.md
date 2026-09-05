# scMarkerAgent (R)

Evidence-grounded cell-type annotation for Seurat `.rds` objects. Every marker claim in
the output resolves to a curated publication and its source sentence.

The companion arm annotates `.h5ad` datasets with the same algorithm. This arm reads
Seurat `.rds` only, and refuses `.h5ad` input rather than converting it.

## Install

From a clone of the repository:

```r
# install.packages("remotes")
remotes::install_local("R", dependencies = TRUE)
```

`dependencies = TRUE` matters: `presto` (the genome-wide differential-expression
engine) is not on CRAN, so it is declared in `Suggests` and resolved through the
`Remotes:` field, which is pinned to a commit of `immunogenomics/presto`. Plain
`install.packages` cannot follow that field; `remotes` or `devtools` can. If `presto`
is missing at run time, `annotate()` stops before the run starts and prints the
one-line install command. The optional report packages (`ggrastr`, `png`, `zip`)
install the same way; without them the pipeline still runs and records what it
skipped.

The pipeline itself -- the same `rflow/` sources the benchmarks ran, with the shared
configuration, prompts, schemas and viewer assets -- ships inside the package under
`inst/scmarkeragent/` and is what `annotate()` runs.

### Pinned environment

The R environment in which the released pipeline was validated is recorded under
`inst/environment/`: `renv.lock` (R 4.5.2; Seurat 5.5.0, SeuratObject 5.4.0,
data.table 1.18.4, Matrix 1.7-4, presto 1.0.0 and the other runtime packages) and
`r-packages.json` (the same versions as a flat list, with restore notes). `DESCRIPTION`
states the looser installation floor (`Seurat >= 5.0.0`); use the lockfile when an
exact reproduction of the released environment is needed
(`renv::restore(lockfile = system.file("environment/renv.lock", package = "scmarkeragent"))`).

## Marker resource

The curated marker resource is required and is distributed separately from the code.

```r
scmarkeragent::download_resources("~/scmarkeragent-resources")
scmarkeragent::verify_resources("~/scmarkeragent-resources")
```

Every file is verified against a shipped SHA-256 index, so a truncated download fails
instead of producing an annotation against a partial database.

## Run

```r
scmarkeragent::annotate(
  input         = "pbmc.rds",
  tag           = "pbmc",
  species       = "Human",
  tissue        = "blood",
  counts_source = "RNA/counts",
  resource_dir  = "~/scmarkeragent-resources"
)
```

| Argument | |
| --- | --- |
| `input` | Seurat `.rds` with raw counts; gene symbols as rownames, unique cell ids |
| `counts_source` | Where the raw counts are, as `<assay>/<layer>`. Required, never inferred |
| `species` | `"Human"`, `"Mouse"` or `"Rat"` |
| `tissue` | Free text, resolved against UBERON |
| `disease` | Default `"Normal"`; separate multiple values with `|` |
| `clustering_resolution` | Leiden resolution. Clustering is always recomputed from counts |
| `resource_dir` | The downloaded marker resource |
| `work_dir` | Cache and results root; defaults to `./.scmarkeragent` |
| `offline` | Run the deterministic stages only, without contacting a model |
| `cluster_annotation` | Set `FALSE` to skip the annotating agent; the label becomes the top of the retrieval order, and every field recording how the label was produced says so |

## Model access

The annotating agent and its quality-control judges call an OpenAI-compatible endpoint.

```r
Sys.setenv(OPENAI_API_KEY = "...")
Sys.setenv(OPENAI_BASE_URL = "...")  # only for a compatible provider
```

Calls are cached by model, reasoning setting and canonical prompt bytes, so rerunning the
same data replays from cache without new requests. Cold calls are written to
`llm_cold_calls.jsonl` with the prompt, its hash, the raw response and the resolved
request settings.

Without credentials the model stages report `skipped_no_credentials` and do not reach the
network. `offline = TRUE` enforces that even when credentials are present.

## Outputs

One directory per run:

| File | Contents |
| --- | --- |
| `cluster_summary.csv` | One row per cluster: the label to read, its parent and subtype, cell state, resolution status, confidence, rationale, quality-control outcome, key markers with detection fractions, the PMCIDs behind them, and any co-occurring identities with their own markers and PMCIDs |
| `marker_evidence.csv` | Every measured marker of every claimed identity: polarity, detection in and out of the cluster, fold change, publication support, PMID/PMCID and source sentence |
| `cluster_evidence.jsonl` | Audit sidecar: every retrieved candidate with its scores and complete measured panel, the agent's turns and tool calls, the exclusion audit, and the judges' verdicts |
| `figures/`, `figure_data/` | Figures and the numbers behind each one |
| `viewer/index.html`, `<tag>_results.zip` | Self-contained interactive report |
| `run_manifest.json`, `resolved_config.json` | Run status, output paths and the effective configuration |

Both arms write the same schema, so results are comparable field by field regardless of
which one produced them.

## Notes

- The pipeline starts from **raw counts** and always recomputes clustering. A partition
  carried by the Seurat object is never reused.
- A cluster carrying more than one population is reported as `mixed`, with the
  co-occurring identities named and their evidence in `cooccurring_markers` /
  `cooccurring_pmcid`. The annotation is not collapsed to a single label.
- One annotation is broadcast to every cell in a cluster, so a minority population inside
  a heterogeneous cluster can be mislabelled.
- A subtype is reported only where the cluster raises at least two markers absent from the
  parent's curated panel. Otherwise the parent is reported and the rejected finer name is
  kept in the evidence table and the claim warnings.
- `annotation_confidence` is an evidence score, not a probability.

## Licence

Noncommercial use only. Code: PolyForm Noncommercial 1.0.0 (`LICENSE`). Marker resource:
CC BY-NC 4.0; its source sentences are quoted from PubMed Central open-access articles
and every row carries its PMCID/PMID. Per-file terms are in the resource's
`resource_manifest.json`.
