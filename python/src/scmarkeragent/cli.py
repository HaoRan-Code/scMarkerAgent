from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bootstrap(argv):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--resource-dir")
    parser.add_argument("--work-dir")
    parser.add_argument("--offline", action="store_true")
    known, _ = parser.parse_known_args(argv)
    if known.resource_dir:
        os.environ["SCMA_RESOURCE_DIR"] = str(Path(known.resource_dir).resolve())
    if known.work_dir:
        os.environ["SCMA_WORK_DIR"] = str(Path(known.work_dir).resolve())
    if known.offline:
        os.environ["SCMA_OFFLINE"] = "1"


def _parser():
    from .config import DEFAULTS

    parser = argparse.ArgumentParser(prog="scmarkeragent")
    sub = parser.add_subparsers(dest="command", required=True)

    annotate = sub.add_parser(
        "annotate",
        help=(
            "Annotate the clusters of one .h5ad dataset against the curated marker "
            "resource: QC, clustering, differential expression, candidate retrieval, "
            "the annotating agent and its judges, and the delivered result package"
        ),
    )
    annotate.add_argument("--input", required=True)
    annotate.add_argument("--tag", default="dataset")
    annotate.add_argument("--species", required=True, choices=["Human", "Mouse", "Rat"])
    annotate.add_argument("--tissue", required=True)
    annotate.add_argument("--disease", default="Normal")
    annotate.add_argument(
        "--counts-source",
        required=True,
        help=(
            "explicit location of the raw counts inside the .h5ad: 'X', 'raw.X', or "
            "'layers/<name>'. Required: the counts matrix is never inferred"
        ),
    )
    annotate.add_argument(
        "--clustering-resolution",
        type=float,
        default=DEFAULTS["preprocessing"]["leiden_resolution"],
        help="Leiden resolution; clustering is always recomputed from the raw counts",
    )
    annotate.add_argument("--resource-dir")
    annotate.add_argument("--work-dir")
    annotate.add_argument("--offline", action="store_true")
    annotate.add_argument(
        "--no-cluster-annotation",
        action="store_true",
        help=(
            "Skip the annotating agent. The run then labels each cluster with the "
            "top of the retrieval order and says so in every field that records how "
            "the label was produced"
        ),
    )

    config = sub.add_parser(
        "config",
        help="Print the effective configuration and the resource paths a run would use",
    )
    config.add_argument("--resource-dir")
    config.add_argument("--work-dir")
    config.add_argument("--offline", action="store_true")

    verify = sub.add_parser(
        "verify-resources",
        help="Check that every resource file a run would read is present, with sizes and hashes",
    )
    verify.add_argument("--resource-dir")
    verify.add_argument("--work-dir")
    verify.add_argument("--offline", action="store_true")

    download = sub.add_parser(
        "download-resources",
        help="Fetch the curated marker resource and verify it against the shipped checksums",
    )
    download.add_argument(
        "--dest",
        required=True,
        help="directory to unpack into; pass it to a run as --resource-dir",
    )
    download.add_argument(
        "--url",
        default="",
        help="archive URL; defaults to the one recorded in the resource index",
    )
    download.add_argument(
        "--check-only",
        action="store_true",
        help="verify an existing directory without downloading anything",
    )
    download.add_argument("--resource-dir")
    download.add_argument("--work-dir")
    download.add_argument("--offline", action="store_true")
    return parser


def _verify_resources():
    from . import config

    paths = {
        "markers": Path(config.DB_CSV),
        "sources": Path(config.DB_SOURCES),
        "uberon_ontology": Path(config.OBO_UBERON),
        "ortholog_human_mouse": Path(config.ORTHO_DIR) / "ortho_Human_to_Mouse.csv",
        "ortholog_human_rat": Path(config.ORTHO_DIR) / "ortho_Human_to_Rat.csv",
    }
    missing = [name for name, path in paths.items() if not path.exists()]
    result = {
        "ok": not missing,
        "missing": missing,
        "files": {
            name: {
                "path": str(path),
                "size": path.stat().st_size if path.exists() else None,
                "sha256": _sha(path) if path.exists() else None,
            }
            for name, path in paths.items()
        },
    }
    print(json.dumps(result, indent=2))
    return 0 if not missing else 2


def _annotate(args):
    if args.offline:
        os.environ["SCMA_OFFLINE"] = "1"

    from . import input_format

    try:
        resolved_input = input_format.check_input(args.input, input_format.PYTHON_ARM)
    except (ValueError, FileNotFoundError) as error:
        raise SystemExit(f"scmarkeragent: {error}") from None

    from . import config
    from . import (
        preprocessing,
        candidate_scoring,
        cluster_annotation,
        llm_client,
        reporting,
    )

    config.ensure_workdirs()
    tag = re.sub(r"[^A-Za-z0-9_.-]+", "_", args.tag)
    output_dir = Path(config.RESULTS_DIR) / tag
    output_dir.mkdir(parents=True, exist_ok=True)
    resolved_path = output_dir / config.DEFAULTS["output"]["resolved_config_file"]
    # One file per run holds every effective setting plus the declared input, which is all a
    # rerun needs. There is deliberately no stage fingerprint or content hash: reuse is decided
    # by the caller, not inferred from a stamp that could silently accept a stale cache.
    config.write_resolved_config(
        resolved_path,
        {
            "tag": tag,
            "input": str(resolved_input),
            "counts_source": args.counts_source,
            "species": args.species,
            "tissue": args.tissue,
            "disease": args.disease,
            "clustering_resolution": float(args.clustering_resolution),
            "cross_species": config.CROSS_SPECIES,
            "compute_umap": config.COMPUTE_UMAP,
            "cluster_annotation": not args.no_cluster_annotation,
            "offline": args.offline,
        },
    )
    manifest = {
        "schema_version": "run-manifest-v1",
        "tag": tag,
        "status": "running",
        "credential_source": llm_client.credential_source(),
        "provider_mode": llm_client.provider_mode(),
        "resolved_config": str(resolved_path),
    }
    manifest_path = output_dir / "run_manifest.json"
    try:
        preprocessing.preprocess(
            tag,
            str(resolved_input),
            args.species,
            args.tissue,
            args.disease.split("|"),
            args.counts_source,
            res=float(args.clustering_resolution),
            cross_species=config.CROSS_SPECIES,
        )
        candidate_scoring.compute_candidate_scores(
            tag,
            top_candidates=int(config.DEFAULTS["retrieval"]["top_candidates"]),
            cross_species=config.CROSS_SPECIES,
        )
        result = cluster_annotation.annotate(
            tag, enabled=not args.no_cluster_annotation
        )
        paths = reporting.generate_report(tag, str(output_dir))
        manifest.update(
            {
                "status": "completed",
                "clusters": len(result),
                "annotation_source": sorted(
                    {record.get("annotation_source", "") for record in result.values()}
                ),
                "llm_status": sorted(
                    {record.get("llm_status", "") for record in result.values()}
                ),
                "outputs": paths,
                "llm_raw_log": str(Path(config.CACHE_DIR) / "llm_cold_calls.jsonl"),
            }
        )
    except Exception as error:
        manifest.update({"status": "failed", "error": repr(error)})
        manifest_path.write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        raise
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    _bootstrap(argv)
    args = _parser().parse_args(argv)
    if args.command == "config":
        from .config import summary

        print(summary())
        return 0
    if args.command == "verify-resources":
        return _verify_resources()
    if args.command == "download-resources":
        from .resources_cli import download_resources, verify_resources

        if args.check_only:
            return verify_resources(args.dest)
        return download_resources(args.dest, args.url)
    return _annotate(args)


if __name__ == "__main__":
    raise SystemExit(main())
