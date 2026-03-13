#!/usr/bin/env python

"""Compute point statistics on spatialdata object for QC."""

import segtraq
import spatialdata as sd
import json
import os
import subprocess

def main():
    print("[START] SegTraQ Point Statistics QC")
    input_path = "${spatialdata_zarr}"
    prefix = "${prefix}"
    centroid_x_key = "${params.segtraq_centroid_x_key}"
    centroid_y_key = "${params.segtraq_centroid_y_key}"
    output_dir = f"segtraq_qc/{prefix}"
    os.makedirs(output_dir, exist_ok=True)
    markers_path = "${markers}"

    #reading the spatial data
    print(f"[INFO] Reading SpatialData object from: {input_path}")
    sdata = sd.read_zarr(input_path)

    #initialiizing segtraq object
    cx_key = centroid_x_key if centroid_x_key not in ("null", "", "None") else None
    cy_key = centroid_y_key if centroid_y_key not in ("null", "", "None") else None
    print("[INFO] Initializing SegTraQ object")

    st = segtraq.SegTraQ(
        sdata,
        tables_centroid_x_key=cx_key,
        tables_centroid_y_key=cy_key
    )

    st.filter_control_and_low_quality_transcripts()

    #reading markers
    with open(markers_path, 'r') as f:
        markers = json.load(f)
        genes_to_test = list(set([g for ct in markers.values() for g in ct.get('positive', []) + ct.get('negative', [])]))

    #computing point statistics
    print(f"[INFO] Computing Point Statistics for QC")
    summary = {}

    if not genes_to_test:
        print("[WARNING] No genes found in markers list. Skipping Point Statistics.")
    else:
        print(f"  Running point statistics for {len(genes_to_test)} genes...")

        st.ps.distance_to_centroid(genes=genes_to_test, restrict_to_within_boundary=True)
        centroid_cols = [f"distance_to_cell_centroid_norm_{g}" for g in genes_to_test
                         if f"distance_to_cell_centroid_norm_{g}" in sdata.tables["table"].obs.columns]
        if centroid_cols:
            mean_dist_cent = float(sdata.tables["table"].obs[centroid_cols].mean().mean())
            summary["mean_normalized_distance_to_centroid"] = mean_dist_cent
            print(f"  Mean normalized distance to centroid: {mean_dist_cent:.4f}")

        st.ps.distance_to_membrane(genes=genes_to_test, restrict_to_within_boundary=True)
        membrane_cols = [f"distance_to_cell_membrane_norm_{g}" for g in genes_to_test
                         if f"distance_to_cell_membrane_norm_{g}" in sdata.tables["table"].obs.columns]
        if membrane_cols:
            mean_dist_memb = float(sdata.tables["table"].obs[membrane_cols].mean().mean())
            summary["mean_normalized_distance_to_membrane"] = mean_dist_memb
            print(f"  Mean normalized distance to membrane: {mean_dist_memb:.4f}")

        if "nucleus_boundaries" in sdata.shapes:
            print("  Calculating compartment localization...")
            st.ps.percentage_transcripts_in_compartments(genes=genes_to_test)

            nuc_cols = [f"pct_nucleus_{g}" for g in genes_to_test if f"pct_nucleus_{g}" in sdata.tables["table"].obs.columns]
            if nuc_cols:
                mean_nuc_pct = float(sdata.tables["table"].obs[nuc_cols].mean().mean())
                summary["mean_percentage_in_nucleus"] = mean_nuc_pct
                print(f"  Mean % in nucleus: {mean_nuc_pct:.2f}%")
        else:
            print("  [SKIP] Nucleus compartment analysis (no nucleus_boundaries found in shapes)")

    #summary
    with open(f"{output_dir}/point_stats_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary written to {output_dir}/point_stats_summary.json")

    obs_csv_path = f"{output_dir}/point_statistics_table.csv"
    sdata.tables["table"].obs.to_csv(obs_csv_path)
    print(f"  Point statistics table saved to {obs_csv_path}")

    version = subprocess.check_output(
        ["pip", "show", "segtraq"], text=True
    )
    segtraq_version = [l for l in version.splitlines() if l.startswith("Version:")][0].split(": ")[1]

    with open("versions.yml", "w") as f:
        f.write(f'"{task.process}":\n')
        f.write(f'  segtraq: "{segtraq_version}"\n')
        f.write(f'  spatialdata: "{sd.__version__}"\n')
    print("[FINISH] SegTraQ Point Statistics QC")

if __name__ == "__main__":
    main()
