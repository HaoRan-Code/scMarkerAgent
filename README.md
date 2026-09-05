# scMarkerAgent

Evidence-grounded cell-type annotation for single-cell RNA-seq, where every marker claim
resolves to a curated publication and its source sentence.

scMarkerAgent does not map cells onto a reference atlas. It retrieves candidate identities
from a curated marker resource, measures each candidate's whole panel in your data, and
asks a language model to decide which identity those measurements support — then puts the
evidence, the competing candidates and the quality-control verdict in the output so a
reader can check the call rather than take it.

## Two arms, one algorithm

The same pipeline is implemented twice so that neither format has to be converted:

| Your data | Package | Reader |
| --- | --- | --- |
| `.h5ad` (AnnData / scanpy) | `python/` | native AnnData |
| `.rds` (Seurat) | `R/` | native Seurat |

Both arms start from **raw counts** and run the same workflow: QC, log-normalization,
HVG/PCA/kNN, Leiden clustering (always recomputed — a partition carried by the input is
never reused), genome-wide Wilcoxon differential expression, candidate retrieval, and the
annotating agent with its quality-control judges. They share one configuration file, one
set of prompts and one marker resource.

An arm refuses the other's format rather than converting it:

```
$ scmarkeragent annotate --input dataset.rds ...
scmarkeragent: .rds input belongs to the R arm, not the PYTHON arm: /data/dataset.rds
  Run it with the R package: scmarkeragent::annotate() instead.
  The two arms read their own object directly; neither converts the other's format.
```

## Install

### Python

```bash
pip install ./python
```

Requires Python 3.12. See `python/pyproject.toml` for pinned dependencies.

### R

```r
# install.packages("remotes")
remotes::install_local("R", dependencies = TRUE)
```

`dependencies = TRUE` matters: `presto` is not on CRAN and is resolved through the
package's `Remotes:` field (pinned to a commit of `immunogenomics/presto`), which
`remotes` or `devtools` can follow and plain `install.packages` cannot. If it is
missing at run time, `annotate()` stops before the run starts and says how to install
it.

## Get the marker resource

The curated marker resource is too large for GitHub and is distributed separately. It is
required — the pipeline has no built-in marker database.

```bash
scmarkeragent download-resources --dest ~/scmarkeragent-resources
```

```r
scmarkeragent::download_resources("~/scmarkeragent-resources")
```

Both verify every file against `required_files/resource_index.json` (SHA-256) after
download. Point a run at the result with `--resource-dir` / `resource_dir =`.

See `required_files/README.md` for the contents, sizes and the archive DOI.

## Quick start

```bash
scmarkeragent annotate \
  --input pbmc.h5ad \
  --tag pbmc \
  --species Human \
  --tissue blood \
  --counts-source X \
  --resource-dir ~/scmarkeragent-resources
```

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

`--counts-source` / `counts_source` is required and never guessed: the pipeline starts
from raw counts, and silently annotating a normalized matrix as if it were counts is the
one failure that produces plausible output from wrong input.

## What a run produces

| File | Contents |
| --- | --- |
| `cluster_summary.csv` | One row per cluster: the label to read, its parent and subtype, resolution status, confidence, rationale, quality-control outcome, key markers with detection fractions, and the PMCIDs behind them |
| `marker_evidence.csv` | Every measured marker of every claimed identity, with detection fractions, fold change, publication support and the source sentence |
| `cluster_evidence.jsonl` | The full audit sidecar: every retrieved candidate with its scores and complete panel, the agent's turns and tool calls, and the judges' verdicts |
| `viewer/index.html` | A self-contained interactive report; also packaged as `<tag>_results.zip` |
| `figures/`, `figure_data/` | Publication figures and the numbers behind each one |

A cluster carrying more than one population is reported as such: `resolution_status` is
`mixed`, the co-occurring identities are named, and `cooccurring_markers` /
`cooccurring_pmcid` carry their evidence — the annotation is not forced into one label.

## Requires a language model

The annotating agent and its quality-control judges call an OpenAI-compatible endpoint.
Set `OPENAI_API_KEY` (and `OPENAI_BASE_URL` for a compatible provider). Every call is
cached by prompt hash, so a rerun of the same data replays without spending anything, and
every cold call is written to an audit log.

`--no-cluster-annotation` runs the retrieval and scoring stages without a model; the label
is then the top of the retrieval order, and every field that records how a label was
produced says so.

## Layout

```
python/          Python package (.h5ad arm)
R/               R package (.rds arm)
required_files/  marker resource — data not in git; see its README for the archive
```

## Licence

Noncommercial use only. Code: PolyForm Noncommercial 1.0.0 (`LICENSE`). Marker resource:
CC BY-NC 4.0; its source sentences are quoted from PubMed Central open-access articles
and every row carries its PMCID/PMID, and the third-party ontology and ortholog files
keep their own terms, recorded per file in the resource's `resource_manifest.json`.
