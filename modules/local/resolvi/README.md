# RESOLVI

Local unsupervised resolVI module using `scviva.model.ResolVI` from scVIVA.
The module runs `bin/resolvi.py` from the pipeline project directory.

## Input

The module accepts `[meta, xenium_dir]`, where `meta.id` identifies the sample
and `xenium_dir` contains:

- `cell_feature_matrix.h5`: Xenium expression counts in 10x HDF5 format.
- `cells.csv.gz`: cell metadata with `cell_id`, `x_centroid` and `y_centroid`.

Coordinates are aligned to expression matrix cell IDs and stored in
`obsm["spatial"]`. Cells with fewer than five total counts are filtered out.
At least 15 cells must remain after filtering.

## Configuration

Pass `--max-epochs` through `task.ext.args` (default: `100`). Training runs on
CPU with `semisupervised=False`. Set `task.ext.prefix` to override the output
prefix, which defaults to `meta.id`.

```groovy
process {
    withName: RESOLVI {
        ext.args = '--max-epochs 100'
        ext.prefix = { "${meta.id}" }
    }
}
```

The module provides a container and a Conda environment with scVIVA 0.1.7.
SpatialData is not required. No packages are installed at task runtime.

## Output

- `h5ad`: `[meta, <prefix>.resolvi.h5ad]`, containing filtered Xenium expression
  data, spatial coordinates and the latent representation in `obsm["X_resolVI"]`.
- `versions`: `versions.yml`, containing the scVIVA version.

The script does not export corrected expression, a saved model or training
history. The input Xenium files are not modified.

## Testing

```bash
nf-test test modules/local/resolvi/tests/main.nf.test --profile docker
python -m unittest discover -s modules/local/resolvi/tests -p 'test_*.py'
```

Use Nextflow >=26.04.6 (the pipeline's minimum) and nf-test >=0.9.5.
The nf-test fixture creates a synthetic Xenium directory with 80 cells and
100 genes and runs one training epoch. Stub execution checks output names and
channels only; its empty H5AD placeholder is not a valid AnnData file.
The Python tests mock scientific dependencies to check CLI arguments, coordinate
alignment, filtering, CPU training and output writing without training a model.