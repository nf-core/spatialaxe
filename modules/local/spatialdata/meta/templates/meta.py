#!/usr/bin/env python

"""Add metadata to SpatialData bundle."""

import spatialdata as sd
import sys
import json
import pandas as pd
import zarr

# Fix zarr v3 + anndata + numcodecs incompatibility:
# anndata's string writer passes numcodecs.VLenUTF8 to zarr.Group.create_array,
# but zarr v3 only accepts ArrayArrayCodec types. OME-Zarr 0.5 requires zarr v3
# for images, so we can't downgrade the store format. Instead, we intercept
# create_array to strip numcodecs codecs and let zarr v3 handle strings natively.
import numcodecs
import zarr.core.group as _zarr_group

_orig_create_array = _zarr_group.Group.create_array


def _v3_compat_create_array(self, *args, **kwargs):
    """Strip numcodecs VLenUTF8 from codec params for zarr v3 compatibility."""
    for param in ("filters", "compressor", "object_codec"):
        val = kwargs.get(param)
        if val is None:
            continue
        if isinstance(val, numcodecs.vlen.VLenUTF8):
            del kwargs[param]
        elif isinstance(val, (list, tuple)):
            cleaned = [v for v in val if not isinstance(v, numcodecs.vlen.VLenUTF8)]
            if len(cleaned) != len(val):
                if cleaned:
                    kwargs[param] = cleaned
                else:
                    del kwargs[param]
    return _orig_create_array(self, *args, **kwargs)


_zarr_group.Group.create_array = _v3_compat_create_array


def _is_arrow_backed(dtype):
    """Check if a pandas dtype is backed by PyArrow."""
    return isinstance(dtype, pd.ArrowDtype) or (
        hasattr(dtype, "storage") and getattr(dtype, "storage", None) == "pyarrow"
    ) or "pyarrow" in str(dtype)


def _convert_df_arrow_to_numpy(df):
    """Convert Arrow-backed dtypes in a DataFrame to numpy object dtype."""
    for col in df.columns:
        dtype = df[col].dtype
        if _is_arrow_backed(dtype):
            df[col] = df[col].astype("object")
        elif isinstance(dtype, pd.CategoricalDtype):
            cats = dtype.categories
            if cats is not None and _is_arrow_backed(cats.dtype):
                df[col] = df[col].cat.rename_categories(cats.astype("object"))
    if _is_arrow_backed(df.index.dtype):
        df.index = pd.Index(df.index.astype("object"))


def convert_arrow_to_numpy(sdata):
    """Convert Arrow-backed dtypes to numpy for anndata zarr write compatibility."""
    for table_key in list(sdata.tables.keys()):
        adata = sdata.tables[table_key]
        _convert_df_arrow_to_numpy(adata.obs)
        _convert_df_arrow_to_numpy(adata.var)


def main():
    print("[START]")

    spatialdata_bundle = "${spatialdata_bundle}"
    xenium_bundle = "${xenium_bundle}"
    output_path = "${prefix}"
    metadata = "${meta}"
    output = "${outputfolder}"

    sdata = sd.read_zarr(f"{spatialdata_bundle}")

    # Convert metadata into dict
    print("[NOTE] Read in provenance ...")
    metadata = metadata.strip("[]")  # Remove square brackets
    pairs = metadata.split(", ")  # Split by comma and space
    metadata = {k: v for k, v in (pair.split(":") for pair in pairs)}  # Create dictionary

    for key in metadata:
        if key not in sdata['raw_table'].uns['spatialdata_attrs']:
            sdata['raw_table'].uns['spatialdata_attrs'][key] = metadata[key]
        else:
            sys.err(f'[ERROR] {key} already exist in sdata[raw_table].uns[spatialdata_attrs].')

    # Add experimental metadata
    print("[NOTE] Read in experiment metadata ...")
    sdata['raw_table'].uns['experiment_xenium'] = ''
    metadata_experiment = f'{xenium_bundle}/experiment.xenium'
    with open(metadata_experiment, "r") as f:
        metadata_experiment = json.load(f)
        sdata['raw_table'].uns['experiment_xenium'] = json.dumps(metadata_experiment)

    # Add gene panel metadata
    print("[NOTE] Read in gene panel metadata ...")
    sdata['raw_table'].uns['gene_panel'] = ''
    metadata_gene_panel = f'{xenium_bundle}/gene_panel.json'
    with open(metadata_gene_panel, "r") as f:
        metadata_gene_panel = json.load(f)
        sdata['raw_table'].uns['gene_panel'] = json.dumps(metadata_gene_panel)

    convert_arrow_to_numpy(sdata)
    sdata.write(f"spatialdata/{output_path}/{output}", overwrite=True, consolidate_metadata=True, sdata_formats=None)

    #Output version information
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'spatialdata: "{sd.__version__}"\\n')

    print("[FINISH]")

if __name__ == "__main__":
    main()
