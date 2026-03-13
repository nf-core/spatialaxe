#!/usr/bin/env python

"""Compute 3D volume metrics on spatialdata object for QC."""

import os
import segtraq
import spatialdata as sd
import json
import subprocess
import pandas as pd


def main():
    print("[START] SegTraQ 3D Volume QC")
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

    #computing metrics
    print(f"[INFO] Computing 3D Volume QC metrics")
    summary = {}

    #heterotypic overlap fraction
    hetero_overlap_df = st.vl.fraction_heterotypic_overlap()
    csv_path = f"{output_dir}/heterotypic_overlap_df.csv"
    hetero_overlap_df.to_csv(csv_path, index=False)
    print(f"  Saved heterotypic overlap stats to {csv_path}")
    mean_area = float(hetero_overlap_df["heterotypic_overlap_area"].mean())
    mean_frac = float(hetero_overlap_df["heterotypic_overlap_fraction"].mean())
    summary["mean_overlap_area"] = mean_area
    summary["mean_overlap_fraction"] = mean_frac
    print(f"  mean_overlap_area: {mean_area}")
    print(f"  mean_overlap_fraction: {mean_frac}")


    #top and bottom z consistency
    sim_top_bottom_df = st.vl.similarity_top_bottom()
    csv_path = f"{output_dir}/similarity_top_bottom.csv"
    sim_top_bottom_df.to_csv(csv_path, index=False)
    print(f"  Saved top-bottom z consistency stats to {csv_path}")
    cosine_sim_top_bottom_z = float(sim_top_bottom_df["cosine_sim_top_bottom_z"].mean())
    summary["cosine_sim_top_bottom_z"] = cosine_sim_top_bottom_z
    print(f"  cosine_sim_top_bottom_z: {cosine_sim_top_bottom_z}")


    #mean VSI per cell
    mean_vsi_df = st.vl.vertical_signal_integrity_per_cell()
    if isinstance(mean_vsi_df, pd.Series):
        mean_vsi_df = mean_vsi_df.to_frame(name="mean_vsi")
    mean_vsi_df.to_csv(f"{output_dir}/mean_vsi.csv", index=False)
    summary["mean_vsi"] = float(mean_vsi_df["mean_vsi"].mean())

    #summary
    with open(f"{output_dir}/volume_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary written to {output_dir}/volume_summary.json")
    print(f"[INFO] Dataframes generated are stored within respective directories")

    version = subprocess.check_output(
        ["pip", "show", "segtraq"], text=True
    )
    segtraq_version = [l for l in version.splitlines() if l.startswith("Version:")][0].split(": ")[1]

    with open("versions.yml", "w") as f:
        f.write(f'"{task.process}":\n')
        f.write(f'  segtraq: "{segtraq_version}"\n')
        f.write(f'  spatialdata: "{sd.__version__}"\n')
    print("[FINISH] SegTraQ 3D Volume QC")

if __name__ == "__main__":
    main()
