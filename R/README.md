# scMarkerAgent (R arm)

Annotates Seurat `.rds` objects. The Python arm in [`../python/`](../python/) annotates
`.h5ad` with the same algorithm; this arm refuses `.h5ad` input rather than converting it.
The step-by-step setup and the smoke test are in the repository [README](../README.md).

## 1. Install

R 4.2 or later. The package is served from its own CRAN-like repository on R-universe
(`https://haoran-code.r-universe.dev`), which also serves `presto`, the genome-wide
differential-expression engine that is not on CRAN, at the commit the release was
validated with. One `install.packages()` therefore brings in everything, as prebuilt
binaries on Windows and macOS:

```r
options(repos = c(scmarkeragent = "https://haoran-code.r-universe.dev",
                  CRAN = "https://cloud.r-project.org"))
install.packages("scmarkeragent")
```

From a clone (after editing the code), set the same `options()` and run
`remotes::install_local("R", dependencies = TRUE)` from the repository root. Seurat,
`presto` and the other runtime packages are declared in `Imports`, so both routes install
them without further flags; `ggrastr`, `png` and `zip` only serve the report and without
them the pipeline still runs and records what it skipped.

Neither route uses the GitHub API, whose anonymous quota (60 calls per hour per IP
address) is what makes `remotes::install_github()` fail with `HTTP error 403` behind
shared addresses; the package therefore carries no `Remotes:` field. If R-universe is
unreachable, `presto` can be installed over the git protocol instead:
`remotes::install_git("https://github.com/immunogenomics/presto.git", ref = "a24772a135c7895a8183b007376050556c60a05b")`.
If `presto` is missing when a run starts, `annotate()` stops before doing anything and
prints both lines. If a CRAN package fails to compile from source, the message names the
missing system library (`libxml2`, `libcurl`, `libpng`, ...); install it with your system's
package manager and rerun.

The pipeline itself (the same `rflow/` sources the benchmarks ran, with the shared
configuration, prompts, schemas and viewer assets) ships under `inst/scmarkeragent/` and
runs in a fresh `Rscript` process. The exact environment the release was validated in is
recorded in `inst/environment/renv.lock` (R 4.5.2, Seurat 5.5.0, presto 1.0.0, ...);
`DESCRIPTION` states the looser floor. For an exact reproduction:
`renv::restore(lockfile = system.file("environment/renv.lock", package = "scmarkeragent"))`.

## 2. Download the marker resource (once)

```r
scmarkeragent::download_resources("~/scmarkeragent-resources")
```

Fetches the archive from Zenodo (about 150 MiB compressed, 760 MiB unpacked), unpacks it
and verifies every file against the shipped SHA-256 index; R's 60-second download timeout
is lifted for the transfer. `verify_resources()` checks an existing directory; `url =`
points at a mirror.

## 3. Set a credential

```r
Sys.setenv(OPENAI_API_KEY = "sk-...")
Sys.setenv(OPENAI_BASE_URL = "https://...")   # only for an OpenAI-compatible provider
Sys.setenv(SCMA_LLM_MODEL = "...")            # optional model override
```

Calls are cached by model, reasoning setting and canonical prompt bytes, so rerunning the
same data replays from cache; cold calls are written to `llm_cold_calls.jsonl`. Without
credentials the model stages report `skipped_no_credentials` and never reach the network;
`offline = TRUE` enforces that even when credentials are present.

## 4. Run

```r
scmarkeragent::annotate(
  input = "pbmc.rds", tag = "pbmc",
  species = "Human", tissue = "blood",
  counts_source = "RNA/counts",
  resource_dir = "~/scmarkeragent-resources"
)
```

| Argument | |
| --- | --- |
| `input` | Seurat `.rds` with raw counts; gene symbols as rownames, unique cell ids |
| `counts_source` | where the raw counts are, as `<assay>/<layer>`. Required, never inferred |
| `species` | `"Human"`, `"Mouse"` or `"Rat"` |
| `tissue` | free text, resolved against UBERON |
| `disease` | default `"Normal"`; several values as a character vector |
| `clustering_resolution` | Leiden resolution, default 0.5. Clustering is always recomputed; a run that yields a single cluster stops with a message, because one-vs-rest differential expression is undefined for it |
| `resource_dir` | the downloaded marker resource, verified before the run starts |
| `work_dir` | cache and results root; default `./.scmarkeragent` |
| `offline` | deterministic stages only, no model contact |
| `cluster_annotation` | `FALSE` skips the annotating agent; the label becomes the top of the retrieval order, and every field recording how the label was produced says so |

`annotate()` returns the parsed `run_manifest.json` invisibly: run status, per-stage
outcomes and the path of every delivered file. `check_input()` tells which arm a file
belongs to.

## 5. Outputs

One directory per run, `<work_dir>/results/<tag>/`:

| File | Contents |
| --- | --- |
| `cluster_summary.csv` | one row per cluster: the label, its parent and subtype, cell state, resolution status, confidence, rationale, QC outcome, key markers with detection fractions, the PMCIDs behind them, and any co-occurring identities with their own markers and PMCIDs |
| `marker_evidence.csv` | every measured marker of every claimed identity: polarity, detection in and out of the cluster, fold change, publication support, PMID/PMCID and source sentence |
| `cluster_evidence.jsonl` | audit sidecar: every retrieved candidate with its scores and complete measured panel, the agent's turns and tool calls, the exclusion audit, the judges' verdicts |
| `figures/`, `figure_data/` | figures and the numbers behind each one |
| `viewer/index.html`, `<tag>_results.zip` | self-contained interactive report |
| `run_manifest.json`, `resolved_config.json` | run status, output paths and the effective configuration |

Both arms write the same schema, so results compare field by field whichever arm produced
them.

## Tests

```r
testthat::test_local("R")   # from the repository root, with the package installed
```

or `R CMD build R && R CMD check scmarkeragent_0.1.0.tar.gz`.

## Licence

Noncommercial use only. Code: PolyForm Noncommercial 1.0.0 (`LICENSE`). Marker resource:
CC BY-NC 4.0; its source sentences are quoted from PubMed Central open-access articles and
every row carries its PMCID/PMID. Per-file terms are in the resource's
`resource_manifest.json`.
