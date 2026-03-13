#!/usr/bin/env python

"""Generate QC plots using SegTraQ plotting module."""

import os
import json
import subprocess
import segtraq
import spatialdata as sd
import scanpy as sc
import matplotlib.pyplot as plt

def main():
    print("[START] SegTraQ Plotting QC")
    input_path = "${spatialdata_zarr}"
    prefix = "${prefix}"
    cell_type_key = "${cell_type_key}"

    centroid_x_key = "${params.segtraq_centroid_x_key}"
    centroid_y_key = "${params.segtraq_centroid_y_key}"

    output_dir = f"segtraq_qc/{prefix}/plots"
    os.makedirs(output_dir, exist_ok=True)

    print(f"[INFO] Reading SpatialData object from: {input_path}")
    sdata = sd.read_zarr(input_path)

    cx_key = centroid_x_key if centroid_x_key not in ("null", "", "None") else None
    cy_key = centroid_y_key if centroid_y_key not in ("null", "", "None") else None
    cell_type_key = cell_type_key if cell_type_key not in ("null", "", "None") else None


    st = segtraq.SegTraQ(
        sdata,
        tables_centroid_x_key=cx_key,
        tables_centroid_y_key=cy_key
    )

    print("[INFO] Preprocessing data for plotting (Normalization, PCA, UMAP)")
    st.filter_control_and_low_quality_transcripts()

    adata = st.sdata.tables[st.tables_key]
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)


    st_dict = {prefix: st}

    print("[INFO] Generating plots...")

    if cell_type_key and cell_type_key in adata.obs.columns:
        print(f"  Plotting cell type proportions using: {cell_type_key}")
        segtraq.pl.celltype_proportions(st_dict, celltype_col=cell_type_key)
        plt.savefig(f"{output_dir}/celltype_proportions.png", bbox_inches='tight')
        plt.close()

        segtraq.pl.umap(st_dict, color=cell_type_key, legend=True)
        plt.savefig(f"{output_dir}/umap_cell_types.png", bbox_inches='tight')
        plt.close()

    if 'transcript_count' in adata.obs.columns:
        segtraq.pl.umap(st_dict, color="transcript_count", legend=True)
        plt.savefig(f"{output_dir}/umap_transcript_count.png", bbox_inches='tight')
        plt.close()

        if cell_type_key and cell_type_key in adata.obs.columns:
            segtraq.pl.boxplot(st_dict, celltype_col=cell_type_key, value_key="transcript_count")
            plt.savefig(f"{output_dir}/boxplot_transcripts_per_type.png", bbox_inches='tight')
            plt.close()

    print(f"[INFO] Plots saved to {output_dir}")

    version = subprocess.check_output(["pip", "show", "segtraq"], text=True)
    segtraq_version = [l for l in version.splitlines() if l.startswith("Version:")][0].split(": ")[1]

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\n')
        f.write(f'  segtraq: "{segtraq_version}"\n')
        f.write(f'  spatialdata: "{sd.__version__}"\n')

    print("[FINISH] SegTraQ Plotting QC")

if __name__ == "__main__":
    main()
