# scMarkerAgent (Python arm)

Annotates `.h5ad` (AnnData) datasets. The R arm in [`../R/`](../R/) annotates Seurat `.rds`
with the same algorithm; this arm refuses `.rds` input rather than converting it. The
step-by-step setup and the smoke test are in the repository [README](../README.md).

## 1. Install

Python 3.12 exactly (`requires-python = ">=3.12,<3.13"`); every dependency is pinned in
`pyproject.toml`.

```bash
conda create -n scmarkeragent python=3.12 -y && conda activate scmarkeragent   # any 3.12 environment works
pip install ./python                     # from a clone of the repository, or:
pip install "git+https://github.com/HaoRan-Code/scMarkerAgent.git#subdirectory=python"
scmarkeragent --help
```

Subcommands: `annotate`, `download-resources`, `verify-resources`, `config`.

## 2. Download the marker resource (once)

```bash
scmarkeragent download-resources --dest ~/scmarkeragent-resources
```

Fetches the archive from Zenodo (about 150 MiB compressed, 760 MiB unpacked), unpacks it
and verifies every file against the shipped SHA-256 index. `--check-only` verifies an
existing directory; `--url` points at a mirror.

## 3. Set a credential

```bash
export OPENAI_API_KEY="sk-..."
export OPENAI_BASE_URL="https://..."   # only for an OpenAI-compatible provider
export SCMA_LLM_MODEL="..."            # optional model override
```

Calls are cached by model, reasoning setting and canonical prompt bytes, so rerunning the
same data replays from cache. Cold calls are written to `llm_cold_calls.jsonl` with the
prompt, its hash, the raw response and the resolved request settings (`SCMA_LLM_RAW_LOG`
moves that file). Without credentials the model stages report `skipped_no_credentials` and
never reach the network; `--offline` enforces that even when credentials are present.

## 4. Run

```bash
scmarkeragent annotate \
  --input pbmc.h5ad --tag pbmc \
  --species Human --tissue blood \
  --counts-source X \
  --resource-dir ~/scmarkeragent-resources
```

| Option | |
| --- | --- |
| `--input` | `.h5ad` with raw counts; gene symbols as `var_names`, unique cell ids |
| `--counts-source` | where the raw counts are: `X`, `raw.X` or `layers/<name>`. Required, never inferred |
| `--species` | `Human`, `Mouse` or `Rat` |
| `--tissue` | free text, resolved against UBERON |
| `--disease` | default `Normal`; several values separated by a vertical bar |
| `--clustering-resolution` | Leiden resolution, default 0.5. Clustering is always recomputed; a run that yields a single cluster stops with a message, because one-vs-rest differential expression is undefined for it |
| `--resource-dir` | the downloaded marker resource |
| `--work-dir` | cache and results root; default `./.scmarkeragent` |
| `--offline` | deterministic stages only, no model contact |
| `--no-cluster-annotation` | skip the annotating agent; the label becomes the top of the retrieval order, and every field recording how the label was produced says so |

`scmarkeragent config` prints the effective configuration; `scmarkeragent verify-resources`
checks the resource files a run would read.

The package also installs `scmarkeragent-r`, which runs the R arm's pipeline (shipped inside
this package under `rflow/`) through `Rscript` when the R dependencies are present. The
documented way to annotate `.rds` data is the R package.

## 5. Outputs

One directory per run, `<work-dir>/results/<tag>/`:

| File | Contents |
| --- | --- |
| `cluster_summary.csv` | one row per cluster: the label, its parent and subtype, cell state, resolution status, confidence, rationale, QC outcome, key markers with detection fractions, the PMCIDs behind them, and any co-occurring identities with their own markers and PMCIDs |
| `marker_evidence.csv` | every measured marker of every claimed identity: polarity, detection in and out of the cluster, fold change, publication support, PMID/PMCID and source sentence |
| `cluster_evidence.jsonl` | audit sidecar: every retrieved candidate with its scores and complete measured panel, the agent's turns and tool calls, the exclusion audit, the judges' verdicts |
| `figures/`, `figure_data/` | figures and the numbers behind each one |
| `viewer/index.html`, `<tag>_results.zip` | self-contained interactive report |
| `run_manifest.json`, `resolved_config.json` | run status, output paths and the effective configuration |

`output_schema.json` in the package is the field-level contract for these files. The output
is reference-free: author labels, purity, ARI and ontology agreement are never emitted.

## Tests

```bash
pip install "./python[test]" && (cd python && pytest)
```

One test verifies a downloaded resource against the shipped manifest; it is skipped unless
`SCMA_RESOURCE_DIR` points at one.

## Licence

Noncommercial use only. Code: PolyForm Noncommercial 1.0.0 (`LICENSE`). Marker resource:
CC BY-NC 4.0; its source sentences are quoted from PubMed Central open-access articles and
every row carries its PMCID/PMID. Per-file terms are in the resource's
`resource_manifest.json`.
