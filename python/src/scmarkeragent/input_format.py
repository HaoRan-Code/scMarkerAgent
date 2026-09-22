
from __future__ import annotations

from pathlib import Path

PYTHON_ARM = "python"
R_ARM = "r"

ARM_BY_SUFFIX = {
    ".h5ad": PYTHON_ARM,
    ".rds": R_ARM,
}

ENTRY_POINT = {
    PYTHON_ARM: "scmarkeragent annotate",
    R_ARM: "the R package: scmarkeragent::annotate()",
}


def arm_for_input(path: str | Path) -> str:
    suffix = Path(path).suffix.lower()
    arm = ARM_BY_SUFFIX.get(suffix)
    if arm is None:
        supported = ", ".join(sorted(ARM_BY_SUFFIX))
        raise ValueError(f"unsupported input format {suffix!r} for {path}; expected one of {supported}")
    return arm


def check_input(path: str | Path, expected_arm: str) -> Path:
    resolved = Path(path).expanduser()
    suffix = resolved.suffix.lower()
    supported = ", ".join(sorted(ARM_BY_SUFFIX))

    if suffix not in ARM_BY_SUFFIX:
        raise ValueError(
            f"unsupported input format {suffix or '(no extension)'!r}: {resolved}\n"
            f"  scMarkerAgent reads {supported} only, and never converts between them.\n"
            f"  Convert the object to one of those formats first, keeping the RAW COUNTS,"
            f" because the pipeline starts from counts."
        )

    arm = ARM_BY_SUFFIX[suffix]
    if arm != expected_arm:
        raise ValueError(
            f"{suffix} input belongs to the {arm.upper()} arm, not the "
            f"{expected_arm.upper()} arm: {resolved}\n"
            f"  Run it with {ENTRY_POINT[arm]} instead.\n"
            f"  The two arms read their own object directly; neither converts the "
            f"other's format."
        )

    if not resolved.exists():
        raise FileNotFoundError(f"input file does not exist: {resolved}")
    if not resolved.is_file():
        raise ValueError(f"input path is not a file: {resolved}")

    return resolved.resolve()
