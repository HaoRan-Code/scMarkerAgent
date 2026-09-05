from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(prog="scmarkeragent-r")
    parser.add_argument("--input", required=True)
    parser.add_argument("--tag", default="dataset")
    parser.add_argument("--species", required=True, choices=["Human", "Mouse", "Rat"])
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--disease", default="Normal")
    parser.add_argument(
        "--counts-source",
        required=True,
        help=(
            "explicit location of the raw counts inside the .rds, as '<assay>/<layer>' "
            "(e.g. 'RNA/counts'). Required: the counts matrix is never inferred"
        ),
    )
    parser.add_argument(
        "--clustering-resolution",
        type=float,
        default=None,
        help="Leiden resolution; clustering is always recomputed from the raw counts. "
        "Defaults to the shared configuration, exactly as the Python CLI does",
    )
    parser.add_argument("--resource-dir", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--offline", action="store_true")
    parser.add_argument(
        "--no-cluster-annotation",
        action="store_true",
        help=(
            "Skip the annotating agent. The run then labels each cluster with the "
            "top of the retrieval order and says so in every field that records how "
            "the label was produced, exactly as the Python CLI does"
        ),
    )
    parser.add_argument("--no-report", action="store_true")
    args = parser.parse_args(argv)

    from . import input_format

    try:
        resolved_input = input_format.check_input(args.input, input_format.R_ARM)
    except (ValueError, FileNotFoundError) as error:
        raise SystemExit(f"scmarkeragent-r: {error}") from None

    rscript = shutil.which("Rscript")
    if not rscript:
        raise SystemExit("Rscript is not available on PATH")
    root = Path(__file__).resolve().parent
    script = root / "rflow" / "run.R"
    env = dict(os.environ)
    env.update(
        {
            "SCMA_INPUT_RDS": str(resolved_input),
            "SCMA_TAG": args.tag,
            "SCMA_SPECIES": args.species,
            "SCMA_TISSUE": args.tissue,
            "SCMA_DISEASE": args.disease,
            "SCMA_COUNTS_SOURCE": args.counts_source,
            "SCMA_RESOURCE_DIR": str(Path(args.resource_dir).resolve()),
            "SCMA_WORK_DIR": str(Path(args.work_dir).resolve()),
            "SCMA_CLUSTER_ANNOTATION": "0" if args.no_cluster_annotation else "1",
            "SCMA_OFFLINE": "1" if args.offline else "0",
            "SCMA_NO_REPORT": "1" if args.no_report else "0",
        }
    )
    if args.clustering_resolution is not None:
        env["SCMA_CLUSTERING_RESOLUTION"] = repr(float(args.clustering_resolution))
    return subprocess.call([rscript, str(script)], env=env)


if __name__ == "__main__":
    raise SystemExit(main())
