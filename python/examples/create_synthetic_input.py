#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite


def main() -> None:
    rng = np.random.default_rng(20260723)
    signature_a = ["CD3D", "CD3E", "TRBC1", "LCK", "IL7R"]
    signature_b = ["MS4A1", "CD79A", "CD37", "CD74", "HLA-DRA"]
    background = [f"SYNTHETIC_GENE_{index:03d}" for index in range(240)]
    genes = signature_a + signature_b + background
    clusters = np.array([0] * 40 + [1] * 40, dtype=np.int32)
    counts = rng.poisson(2.5, size=(len(clusters), len(genes))).astype(np.int32)
    counts[:40, : len(signature_a)] += rng.poisson(8, size=(40, len(signature_a)))
    counts[40:, len(signature_a) : len(signature_a) + len(signature_b)] += rng.poisson(
        8, size=(40, len(signature_b))
    )
    cells = [f"synthetic_cell_{index:03d}" for index in range(len(clusters))]
    obs = pd.DataFrame(
        {"cluster": clusters},
        index=pd.Index(np.asarray(cells, dtype=object), dtype=object),
    )
    var = pd.DataFrame(index=pd.Index(np.asarray(genes, dtype=object), dtype=object))
    data = ad.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)
    output = Path(__file__).with_name("synthetic_input.h5ad")
    data.write_h5ad(output)
    root = output.parent
    mmwrite(root / "synthetic_counts.mtx", sp.csr_matrix(counts).T)
    (root / "synthetic_genes.tsv").write_text("\n".join(genes) + "\n", encoding="utf-8")
    (root / "synthetic_cells.tsv").write_text("\n".join(cells) + "\n", encoding="utf-8")
    pd.DataFrame({"cell": cells, "cluster": clusters}).to_csv(
        root / "synthetic_clusters.tsv", sep="\t", index=False
    )
    print(output)


if __name__ == "__main__":
    main()
