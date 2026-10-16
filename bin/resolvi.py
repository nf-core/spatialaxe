
#!/usr/bin/env python3

import argparse
from pathlib import Path

import scanpy as sc
import pandas as pd
from scviva.model import ResolVI

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--max-epochs", type=int, default=100)
args = parser.parse_args()

# Load Xenium expression matrix
adata = sc.read_10x_h5(args.input / "cell_feature_matrix.h5")

# Load cell coordinates
cells = pd.read_csv(
    args.input / "cells.csv.gz",
    usecols=["cell_id", "x_centroid", "y_centroid"]
).set_index("cell_id")

adata.obsm["spatial"] = cells.loc[
    adata.obs_names, ["x_centroid", "y_centroid"]
].to_numpy()

print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

# Filter low-count cells
sc.pp.filter_cells(adata, min_counts=5)

if adata.n_obs < 15:
    raise ValueError("Too few cells for resolVI.")

# Prepare and train resolVI
ResolVI.setup_anndata(
    adata,
    prepare_data_kwargs={"spatial_rep": "spatial"}
)

model = ResolVI(adata, semisupervised=False)
model.train(max_epochs=args.max_epochs, accelerator="cpu")

# Store latent representation
adata.obsm["X_resolVI"] = model.get_latent_representation()

# Save output
adata.write_h5ad(args.output)
print(f"SUCCESS: {args.output}")
