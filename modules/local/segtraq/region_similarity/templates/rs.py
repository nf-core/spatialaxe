#!/usr/bin/env python

"""Compute region similarity metrics on spatialdata object for QC."""

import os
import segtraq
import spatialdata as sd
import json
import subprocess
import pandas as pd


def main():
    print("[START] SegTraQ Region Similarity QC")
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
        points_background_id = 0,
        tables_centroid_x_key= cx_key,
        tables_centroid_y_key= cy_key,
    )

    print(f"[INFO] checking the presence of cell and nuclear masks")
    if hasattr(st.sdata, "shapes") and st.sdata.shapes:
        print("Found shapes. Moving ahead.")
    else:
        raise ValueError("A nuclear and cell segmentation object in 'shapes' is required.")

    #computing metrics
    print(f"[INFO] Computing region similarity QC metrics")
    summary = {}

    #match nuclei to cells
    nuclei_to_cells_df = st.rs.match_nuclei_to_cells()
    csv_path = f"{output_dir}/match_nucleus_to_cell.csv"
    nuclei_to_cells_df.to_csv(csv_path, index=False)
    print(f"  Saved full nuclei matching stats to {csv_path}")
    mean_iou = float(nuclei_to_cells_df["iou"].mean())
    summary["mean_nucleus_cell_iou"] = mean_iou
    print(f"  mean nucleus-cell IoU: {mean_iou}")

    #similarity between border and neighbourhood
    sim_border_neighborhood_df = st.rs.similarity_border_neighborhood()
    csv_path = f"{output_dir}/similarity_border_neighbourhood.csv"
    sim_border_neighborhood_df.to_csv(csv_path, index=False)
    print(f"  Saved similarity between border and neighborhood stats to {csv_path}")
    mean_sim_cent_border = float(sim_border_neighborhood_df["similarity_center_border"].mean())
    summary["mean_similarity_center_border"] = mean_sim_cent_border
    print(f"  mean_similarity_center_border: {mean_sim_cent_border}")
    mean_sim_border_neigh = float(sim_border_neighborhood_df["similarity_border_neighborhood"].mean())
    summary["mean_similarity_border_neighborhood"] = mean_sim_border_neigh
    print(f"  mean_similarity_border_neighborhood: {mean_sim_border_neigh}")
    mean_ratio = float(sim_border_neighborhood_df["ratio_border_neighborhood_to_center"].mean())
    summary["mean_ratio_border_neighborhood_to_center"] = mean_ratio
    print(f"  mean_ratio_border_neighborhood_to_center: {mean_ratio}")

    #similarity between nucleus and cell
    sim_nucleus_cell = st.rs.similarity_nucleus_cell()
    if isinstance(sim_nucleus_cell, pd.Series):
        sim_nucleus_cell = sim_nucleus_cell.to_frame(name="similarity_nucleus_cell")
    sim_nucleus_cell.to_csv(f"{output_dir}/similarity_nucleus_cell.csv", index=False)
    summary["mean_similarity_nucleus_cell"] = float(sim_nucleus_cell.iloc[:, -1].mean())

    #similarity between nucleus and cytoplasm
    sim_nucleus_cyto = st.rs.similarity_nucleus_cytoplasm()
    if isinstance(sim_nucleus_cyto, pd.Series):
        sim_nucleus_cyto = sim_nucleus_cyto.to_frame(name="similarity_nucleus_cytoplasm")
    sim_nucleus_cyto.to_csv(f"{output_dir}/similarity_nucleus_cytoplasm.csv", index=False)
    summary["mean_similarity_nucleus_cytoplasm"] = float(sim_nucleus_cyto.iloc[:, -1].mean())

    #summary
    with open(f"{output_dir}/region_similarity_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary written to {output_dir}/region_similarity_summary.json")
    print(f"[INFO] Dataframes generated are stored within respective directories")

    version = subprocess.check_output(
        ["pip", "show", "segtraq"], text=True
    )
    segtraq_version = [l for l in version.splitlines() if l.startswith("Version:")][0].split(": ")[1]

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\n')
        f.write(f'  segtraq: "{segtraq_version}"\n')
        f.write(f'  spatialdata: "{sd.__version__}"\n')
    print("[FINISH] SegTraQ Region Similarity QC")

if __name__ == "__main__":
    main()
