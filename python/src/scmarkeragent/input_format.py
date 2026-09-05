"""Which execution arm reads a given input format.

The two arms are not interchangeable back-ends for one reader: each reads its native
object directly, so an `.h5ad` dataset is annotated by the Python arm and a Seurat
`.rds` dataset by the R arm, and neither format is converted to the other first. That
rule is stated once here so the entry points and the batch tooling cannot drift apart
on it.
"""

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
    """The arm that reads this input natively, by file extension.

    Raises for any other extension rather than guessing: silently routing an
    unrecognized format to one arm would either fail deep inside a reader or, worse,
    succeed against a partially parsed object.
    """
    suffix = Path(path).suffix.lower()
    arm = ARM_BY_SUFFIX.get(suffix)
    if arm is None:
        supported = ", ".join(sorted(ARM_BY_SUFFIX))
        raise ValueError(f"unsupported input format {suffix!r} for {path}; expected one of {supported}")
    return arm


def check_input(path: str | Path, expected_arm: str) -> Path:
    """Validate one entry point's input before any stage runs, or explain why not.

    The check is at the entry point rather than inside the reader because the two
    failures it prevents look nothing alike to a user. A `.rds` handed to the Python arm
    dies inside an HDF5 parser with a message about a corrupt file, and a path typo dies
    even later, after preprocessing has already written a cache directory. Both are the
    same mistake -- the wrong file, or the wrong arm for the right file -- and both are
    answerable here, in one sentence, before anything is computed.
    """
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
