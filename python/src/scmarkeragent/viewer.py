"""Build a portable, reference-free viewer for one production result.

The viewer's CSS, JS and HTML skeleton live in `resources/static/viewer/`, SHARED with
the R arm's `rflow/viewer.R`: both arms splice the same bytes around the same payload,
which is what makes the delivered `index.html` identical no matter which arm wrote it.
"""

from __future__ import annotations

import csv
import html
import json
import math
import shutil
import zipfile
from datetime import UTC, datetime
from pathlib import Path

_STATIC = Path(__file__).resolve().parent / "resources" / "static" / "viewer"
CSS = (_STATIC / "viewer.css").read_text(encoding="utf-8")
JS = (_STATIC / "viewer.js").read_text(encoding="utf-8")
_TEMPLATE = (_STATIC / "index_template.html").read_text(encoding="utf-8")

def _read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_jsonl(path: Path) -> list[dict]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line
    ]


def _safe_json(value) -> str:
    def clean(item):
        if isinstance(item, dict):
            return {str(key): clean(entry) for key, entry in item.items()}
        if isinstance(item, list):
            return [clean(entry) for entry in item]
        if isinstance(item, float) and (math.isnan(item) or math.isinf(item)):
            return None
        return item

    return (
        json.dumps(
            clean(value), ensure_ascii=False, separators=(",", ":"), allow_nan=False
        )
        .replace("</", "<\\/")
        .replace("<!--", "<\\!--")
    )


def _copy_viewer_files(result_dir: Path, viewer_dir: Path) -> list[str]:
    assets = viewer_dir / "assets" / "figures"
    downloads = viewer_dir / "downloads"
    assets.mkdir(parents=True)
    downloads.mkdir(parents=True)
    copied = []
    figure_names = (
        "dotplot_celltype_markers.png",
        "dotplot_celltype_markers_by_celltype.png",
        "dotplot_celltype_markers_publication_support.png",
        "umap_celltype_labeled.png",
        "umap_cluster.png",
    )
    for name in figure_names:
        source = result_dir / "figures" / name
        if source.is_file():
            target = assets / name
            shutil.copy2(source, target)
            copied.append(str(target.relative_to(viewer_dir)))
    for name in (
        "cluster_summary.csv",
        "marker_evidence.csv",
        "cluster_evidence.jsonl",
        "markers_all_by_cluster.csv",
        "markers_significant_by_cluster.csv",
    ):
        source = result_dir / name
        if not source.is_file():
            continue
        target = downloads / name
        shutil.copy2(source, target)
        copied.append(str(target.relative_to(viewer_dir)))
    for source in sorted((result_dir / "figures").glob("*.pdf")):
        target = downloads / source.name
        shutil.copy2(source, target)
        copied.append(str(target.relative_to(viewer_dir)))
    return copied


def _figure_html(viewer_dir: Path) -> tuple[str, str]:
    figures = [
        ("identity-markers", "Identity markers", "dotplot_celltype_markers.png"),
        (
            "identity-markers-celltype",
            "Identity markers by cell type",
            "dotplot_celltype_markers_by_celltype.png",
        ),
        (
            "identity-markers-support",
            "Publication support",
            "dotplot_celltype_markers_publication_support.png",
        ),
        ("cell-types", "Cell types", "umap_celltype_labeled.png"),
        ("clusters", "Clusters", "umap_cluster.png"),
    ]
    available = [
        item
        for item in figures
        if (viewer_dir / "assets" / "figures" / item[2]).is_file()
    ]
    tabs = "".join(
        f'<button class="tab figure-tab{" active" if index == 0 else ""}" '
        f'data-figure="{identifier}">{html.escape(label)}</button>'
        for index, (identifier, label, _) in enumerate(available)
    )
    images = "".join(
        f'<img id="{identifier}" class="{"active" if index == 0 else ""}" '
        f'src="assets/figures/{html.escape(filename)}" alt="{html.escape(label)}">'
        for index, (identifier, label, filename) in enumerate(available)
    )
    return tabs, images


def _build_html(tag: str, result_dir: Path, viewer_dir: Path) -> str:
    clusters = _read_csv(result_dir / "cluster_summary.csv")
    markers = _read_csv(result_dir / "marker_evidence.csv")
    evidence = _read_jsonl(result_dir / "cluster_evidence.jsonl")
    total_cells = sum(int(row["n_cells"]) for row in clusters)
    resolutions = sorted(
        {row["resolution_status"] for row in clusters if row["resolution_status"]}
    )
    options = '<option value="">All resolution states</option>' + "".join(
        f'<option value="{html.escape(value)}">{html.escape(value)}</option>'
        for value in resolutions
    )
    figure_tabs, figure_images = _figure_html(viewer_dir)
    data = _safe_json({"clusters": clusters, "markers": markers, "evidence": evidence})
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(tag)} — scMarkerAgent</title>
<style>{CSS}</style>
</head>
<body>
<header class="topbar"><div class="topbar-inner">
  <div class="brand">scMarkerAgent</div>
  <div class="contract">Reference-free production result · one dataset</div>
</div></header>
<main>
  <h1>{html.escape(tag)}</h1>
  <p class="lead">Marker-grounded cluster annotations with measured, LLM-validated,
  source-backed cell-type identity markers.</p>
  <div class="notice"><strong>Reference-free production output.</strong>
  No author annotation, purity, reference metric, or secondary label is read or displayed.</div>
  <div class="cards">
    <div class="card"><div class="metric-label">Cells</div><div class="metric-value">{total_cells:,}</div><div class="metric-note">after production QC</div></div>
    <div class="card"><div class="metric-label">Clusters</div><div class="metric-value">{len(clusters)}</div><div class="metric-note">one reported annotation per cluster</div></div>
    <div class="card"><div class="metric-label">Identity markers</div><div class="metric-value">{len(markers):,}</div><div class="metric-note">LLM-validated and DB-cited</div></div>
  </div>
  <div class="downloads">
    <a class="button" href="downloads/cluster_summary.csv">Cluster summary</a>
    <a class="button" href="downloads/marker_evidence.csv">Identity marker evidence</a>
    <a class="button" href="downloads/markers_all_by_cluster.csv">All markers per cluster</a>
    <a class="button" href="downloads/markers_significant_by_cluster.csv">Significant markers per cluster</a>
    <a class="button" href="downloads/cluster_evidence.jsonl">Structured evidence</a>
  </div>
  <section class="panel"><h2>Figures</h2><div class="tabs">{figure_tabs}</div>
    <div class="figure-stage">{figure_images}</div>
  </section>
  <section class="panel"><h2>Cluster explorer</h2>
    <div class="toolbar">
      <input id="cluster-search" type="search" placeholder="Search cluster, annotation, or marker">
      <select id="resolution-filter">{options}</select>
      <span id="cluster-count" class="count"></span>
    </div>
    <div class="table-wrap"><table><thead><tr>
      <th class="num"><button data-sort="cluster_id">Cluster</button></th>
      <th class="num"><button data-sort="n_cells">Cells</button></th>
      <th><button data-sort="cell_type_annotation">Cell type</button></th>
      <th><button data-sort="cell_subtype_annotation">Cell subtype</button></th>
      <th><button data-sort="cell_ontology">Cell Ontology</button></th>
      <th><button data-sort="annotation_confidence">Confidence</button></th>
      <th><button data-sort="resolution_status">Resolution</button></th>
      <th><button data-sort="annotation_source">Decided by</button></th>
    </tr></thead><tbody id="cluster-body"></tbody></table></div>
  </section>
  <section id="cluster-detail" class="panel detail"><div class="empty">Select a cluster.</div></section>
  <p class="footer">Confidence is the high / medium / low judgement the annotating agent made
  from the supplied evidence. Every retrieved candidate, its independent review, and the three
  relative scores behind the retrieval order are kept in the structured evidence file. The
  Cell Ontology column is displayed as curated and is never used to decide, merge or rank an
  annotation; N/A means the curated free-text cell type carries no ontology id.</p>
</main>
<script id="result-data" type="application/json">{data}</script>
<script>{JS}</script>
</body>
</html>
"""


def _write_zip(viewer_dir: Path, output: Path) -> None:
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(viewer_dir.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(viewer_dir).as_posix())


def build_dataset_viewer(tag: str, outdir: str | Path) -> dict[str, str]:
    result_dir = Path(outdir)
    viewer_dir = result_dir / "viewer"
    if viewer_dir.exists():
        shutil.rmtree(viewer_dir)
    viewer_dir.mkdir(parents=True)
    copied = _copy_viewer_files(result_dir, viewer_dir)
    html_path = viewer_dir / "index.html"
    html_path.write_text(
        _build_html(tag, result_dir, viewer_dir),
        encoding="utf-8",
    )
    copied.append("index.html")
    manifest = {
        "schema_version": "scmarkeragent-dataset-viewer-v1",
        "dataset": tag,
        "reference_free": True,
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "entrypoint": "index.html",
        "files": sorted(copied),
    }
    manifest_path = viewer_dir / "viewer_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    zip_path = result_dir / f"{tag}_results.zip"
    if zip_path.exists():
        zip_path.unlink()
    _write_zip(viewer_dir, zip_path)
    return {
        "viewer": str(html_path.resolve()),
        "viewer_zip": str(zip_path.resolve()),
    }
