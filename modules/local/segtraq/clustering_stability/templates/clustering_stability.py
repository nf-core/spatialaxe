#!/usr/bin/env python

"""Compute clustering stability metrics on spatialdata object for QC."""

import os
import segtraq
import spatialdata as sd
import json
import subprocess
import scanpy as sc


def main():
    print("[START] SegTraQ Clustering Stability QC")
    input_path = "${spatialdata_zarr}"
    prefix = "${prefix}"
    centroid_x_key = "${params.segtraq_centroid_x_key}"
    centroid_y_key = "${params.segtraq_centroid_y_key}"
    output_dir = f"segtraq_qc/{prefix}"
    os.makedirs(output_dir, exist_ok=True)

    #reading the spatial data
    print(f"[INFO] Reading SpatialData object from: {input_path}")
    sdata = sd.read_zarr(input_path)

    #initialiizing segtraq object
    cx_key = centroid_x_key if centroid_x_key not in ("null", "", "None") else None
    cy_key = centroid_y_key if centroid_y_key not in ("null", "", "None") else None
    print("[INFO] Initializing SegTraQ object")
    st = segtraq.SegTraQ(
        sdata,
        images_key = None,
        tables_area_key = None,
        points_background_id =0,
        tables_centroid_x_key=cx_key,
        tables_centroid_y_key=cy_key,
    )
    #normalizing and log-transforming for clustering stability metrics
    print(f"[INFO] Normalizing data for clustering stability metrics")
    adata = st.sdata.tables["table"]
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    # normalizing and log-transforming the counts
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    # computing a PCA and neighbors
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)


    #computing metrics
    print(f"[INFO] Computing clustering stability QC metrics")
    summary = {}

    #adjusted_rand_index
    adjusted_rand_index = st.cs.adjusted_rand_index()
    summary["adjusted_rand_index"] = float(adjusted_rand_index)
    print(f"  adjusted_rand_index: {adjusted_rand_index}")

    #cluster_connectedness
    cluster_connectedness = st.cs.cluster_connectedness(use_weights=True)
    summary["cluster_connectedness"] = float(cluster_connectedness)
    print(f" cluster_connectedness: {cluster_connectedness}")

    #purity
    purity = st.cs.purity()
    summary["purity"] = float(purity)
    print(f" purity: {purity}")

    #silhouette_score
    silhouette_score = st.cs.silhouette_score()
    summary["silhouette_score"] = float(silhouette_score)
    print(f" silhouette_score: {silhouette_score}")

    #summary
    with open(f"{output_dir}/clustering_stability_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary written to {output_dir}/clustering_stability_summary.json")

    version = subprocess.check_output(
        ["pip", "show", "segtraq"], text=True
    )
    segtraq_version = [l for l in version.splitlines() if l.startswith("Version:")][0].split(": ")[1]

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\n')
        f.write(f'  segtraq: "{segtraq_version}"\n')
        f.write(f'  spatialdata: "{sd.__version__}"\n')
    print("[FINISH] SegTraQ Clustering Stability QC")

if __name__ == "__main__":
    main()
