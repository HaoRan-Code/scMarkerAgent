# scMarkerAgent (Python)

Evidence-grounded cell-type annotation for `.h5ad` datasets. Every marker claim in the
output resolves to a curated publication and its source sentence.

The companion R package annotates Seurat `.rds` objects with the same algorithm. This arm
refuses `.rds` input rather than converting it.

## Install

From a clone of the repository:

```bash
pip install ./python
```

Python 3.12. Dependencies are pinned in `pyproject.toml`.

## Marker resource

The curated marker resource is required and is distributed separately from the code.

```bash
scmarkeragent download-resources --dest ~/scmarkeragent-resources
scmarkeragent download-resources --check-only --dest ~/scmarkeragent-resources
```

Every file is verified against a shipped SHA-256 index, so a truncated download fails
instead of producing an annotation against a partial database.

## Run

```bash
scmarkeragent annotate \
  --input pbmc.h5ad \
  --tag pbmc \
  --species Human \
  --tissue blood \
  --counts-source X \
  --resource-dir ~/scmarkeragent-resources
```

| Option | |
| --- | --- |
| `--input` | `.h5ad` file with raw counts; gene symbols as `var_names`, unique cell ids |
| `--counts-source` | Where the raw counts are: `X`, `raw.X`, or `layers/<name>`. Required, never inferred |
| `--species` | `Human`, `Mouse` or `Rat` |
| `--tissue` | Free text, resolved against UBERON |
| `--disease` | Default `Normal`; separate multiple values with `|` |
| `--clustering-resolution` | Leiden resolution. Clustering is always recomputed from counts |
| `--resource-dir` | The downloaded marker resource |
| `--work-dir` | Cache and results root; defaults to `./.scmarkeragent` |
| `--offline` | Run the deterministic stages only, without contacting a model |
| `--no-cluster-annotation` | Skip the annotating agent; the label becomes the top of the retrieval order, and every field recording how the label was produced says so |

Other subcommands: `config` prints the effective configuration, `verify-resources` checks
the resources a run would use.

## Model access

The annotating agent and its quality-control judges call an OpenAI-compatible endpoint.

```bash
export OPENAI_API_KEY="..."
export OPENAI_BASE_URL="..."   # only for a compatible provider
export SCMA_LLM_MODEL="..."    # optional model override
```

Calls are cached by model, reasoning setting and canonical prompt bytes, so rerunning the
same data replays from cache without new requests. Cold calls are written to
`llm_cold_calls.jsonl` with the prompt, its hash, the raw response and the resolved
request settings; `SCMA_LLM_RAW_LOG` moves that file.

Without credentials the model stages report `skipped_no_credentials` and do not reach the
network. `--offline` enforces that even when credentials are present.

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

`output_schema.json` in the package is the field-level contract for these files.

## Notes

- The pipeline starts from **raw counts** and always recomputes clustering. A partition
  carried by the input object is never reused.
- A cluster carrying more than one population is reported as `mixed`, with the
  co-occurring identities named and their evidence in `cooccurring_markers` /
  `cooccurring_pmcid`. The annotation is not collapsed to a single label.
- One annotation is broadcast to every cell in a cluster, so a minority population inside
  a heterogeneous cluster can be mislabelled.
- A subtype is reported only where the cluster raises at least two markers absent from the
  parent's curated panel. Otherwise the parent is reported and the rejected finer name is
  kept in the evidence table and the claim warnings.
- `annotation_confidence` is an evidence score, not a probability.
- The output is reference-free: author labels, purity, ARI and ontology agreement are
  never emitted.

## Licence

Noncommercial use only. Code: PolyForm Noncommercial 1.0.0 (`LICENSE`). Marker resource:
CC BY-NC 4.0; its source sentences are quoted from PubMed Central open-access articles
and every row carries its PMCID/PMID. Per-file terms are in the resource's
`resource_manifest.json`.
