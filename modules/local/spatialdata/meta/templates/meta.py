#!/usr/bin/env python

"""Add metadata to SpatialData bundle."""

import spatialdata as sd
import sys
import json

def main():
    print("[START]")

    spatialdata_bundle = "${spatialdata_bundle}"
    xenium_bundle = "${xenium_bundle}"
    metadata = "${meta}"
    output = "spatialdata_meta"

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

    sdata.write(f"./{output}", overwrite=True, consolidate_metadata=True, format=None)

    #Output version information
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'spatialdata: "{sd.__version__}"\\n')

    print("[FINISH]")

if __name__ == "__main__":
    main()
