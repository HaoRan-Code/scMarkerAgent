# Marker resource

The curated marker resource scMarkerAgent annotates against. It is **required** — neither
package ships a marker database — and it is too large for git, so it lives in the archive
below; a fresh clone carries only this README, the checksum index (`resource_index.json`,
which both packages ship a copy of) and the bundle manifest.

## Download

```bash
scmarkeragent download-resources --dest ~/scmarkeragent-resources
```

```r
scmarkeragent::download_resources("~/scmarkeragent-resources")
```

Both read `resource_index.json` in this directory and verify every file's SHA-256 after
download. A partial or corrupted download therefore fails loudly instead of producing an
annotation against half a database.

To install it by hand, download the archive, unpack it anywhere, and pass that directory
as `--resource-dir` (Python) or `resource_dir =` (R).

## Contents

| File | Size | What it is |
| --- | ---: | --- |
| `scmarkeragent_curated_sources.csv` | 612 MiB | The source sentence behind every marker record, with PMID/PMCID — row-for-row 1:1 with the marker table. This is what makes a call checkable |
| `scmarkeragent_curated_markers.csv` | 136 MiB | Curated marker–cell-type pairs with polarity, publication support and evidence tier |
| `ontology/uberon-basic.obo` | ~11 MiB | Tissue context resolution |
| `ortholog/ortho_Human_to_*.csv` | 419 KiB | Cross-species marker pooling |
| `resource_manifest.json` | 2 KiB | Size, checksum and licence of every file above |

Six files, ~760 MiB, all of them verified by `resource_index.json`.

## Archive

The bundle is archived on Zenodo as `scmarkeragent-curated.tar.gz` under the concept DOI
[doi:10.5281/zenodo.22333283](https://doi.org/10.5281/zenodo.22333283), which always resolves
to the latest version of the record. Both commands download from the version-pinned file URL
recorded as `archive_url` in `resource_index.json` (the version whose size and SHA-256 the
index carries), so the bytes verified on arrival are exactly the bytes the index describes.
To fetch the bundle from a different location, pass the URL explicitly:

```bash
scmarkeragent download-resources --dest ~/res --url https://example.org/bundle.tar.gz
```

```r
scmarkeragent::download_resources("~/res", url = "https://example.org/bundle.tar.gz")
```

## Licence

Noncommercial use only (CC BY-NC 4.0). The source sentences are quoted from PubMed
Central open-access articles and every row carries its PMCID/PMID. The ontology and
ortholog files keep their own third-party terms, recorded per file in
`resource_manifest.json`.
