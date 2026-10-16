#!/usr/bin/env python3

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def make_xenium(directory="synthetic"):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)
    counts = sparse.csc_matrix(rng.poisson(np.linspace(1, 12, 100), size=(80, 100)).T)
    cell_ids = [f"cell_{index}" for index in range(80)]
    gene_ids = [f"gene_{index}" for index in range(100)]
    with h5py.File(directory / "cell_feature_matrix.h5", "w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset("data", data=counts.data)
        matrix.create_dataset("indices", data=counts.indices)
        matrix.create_dataset("indptr", data=counts.indptr)
        matrix.create_dataset("shape", data=counts.shape)
        matrix.create_dataset("barcodes", data=np.asarray(cell_ids, dtype="S"))
        features = matrix.create_group("features")
        features.create_dataset("id", data=np.asarray(gene_ids, dtype="S"))
        features.create_dataset("name", data=np.asarray(gene_ids, dtype="S"))
        features.create_dataset("feature_type", data=np.asarray(["Gene Expression"] * 100, dtype="S"))
        features.create_dataset("genome", data=np.asarray(["synthetic"] * 100, dtype="S"))
    pd.DataFrame(
        {
            "cell_id": cell_ids,
            "x_centroid": rng.uniform(0, 100, 80),
            "y_centroid": rng.uniform(0, 100, 80),
        }
    ).iloc[::-1].to_csv(directory / "cells.csv.gz", index=False)


if __name__ == "__main__":
    make_xenium()