#!/usr/bin/env python

"""Compute supervised metrics on spatialdata object for QC."""

import os
import segtraq
import spatialdata as sd
import json
import subprocess
import scanpy as sc
import pandas as pd


def main():
    print("[START] SegTraQ Supervised metrics QC")
    input_path = "${spatialdata_zarr}"
    prefix = "${prefix}"
    cell_type_key = "${cell_type_key}"
    markers_path = "${markers}"

    centroid_x_key = "${params.segtraq_centroid_x_key}"
    centroid_y_key = "${params.segtraq_centroid_y_key}"
    output_dir = f"segtraq_qc/{prefix}"
    os.makedirs(output_dir, exist_ok=True)

    #loading markers
    with open(markers_path, 'r') as f:
        markers = json.load(f)

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

    #normalization
    adata = sdata.tables[st.tables_key]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)


    #computing metrics
    print(f"[INFO] Computing supervised QC metrics")
    summary = {}

    #marker purity
    marker_purity_df = st.sp.marker_purity(
        cell_type_key=cell_type_key,
        markers=markers
    )
    csv_path = f"{output_dir}/marker_purity.csv"
    marker_purity_df.to_csv(csv_path, index=False)
    print(f"  Saved marker purity stats to {csv_path}")
    mean_positive_F1 = float(marker_purity_df["positive_F1"].mean())
    summary["mean_positive_F1"] = mean_positive_F1
    print(f"  mean_positive_F1: {mean_positive_F1}")
    mean_negative_F1 = float(marker_purity_df["negative_F1"].mean())
    summary["mean_negative_F1"] = mean_negative_F1
    print(f"  mean_negative_F1: {mean_negative_F1}")
    mean_F1_purity = float(marker_purity_df["F1_purity"].mean())
    summary["F1_purity"] = mean_F1_purity
    print(f"  mean_F1_purity: {mean_F1_purity}")

    #mutually_exclusive_coexpression_rate
    mecr_df = st.sp.mutually_exclusive_coexpression_rate(
            markers=markers
        )
    csv_path = f"{output_dir}/mutually_exclusive_coexpression_rate.csv"
    mecr_df.to_csv(csv_path, index=False)
    print(f"  Saved mutually exclusive coexpression rate stats to {csv_path}")
    sig_count = (mecr_df["pvalue"] < 0.05).sum()
    total_pairs = len(mecr_df)
    summary["total_gene_pairs_tested"] = int(total_pairs)
    summary["significant_exclusive_pairs_count"] = int(sig_count)
    summary["percentage_significant_exclusivity"] = float((sig_count / total_pairs) * 100) if total_pairs > 0 else 0.0

    #neighbor contamination
    per_cell_df, strength_df, binary_df = st.sp.neighbor_contamination(
        cell_type_key=cell_type_key,
        markers=markers
    )
    csv_path = f"{output_dir}/per_cell_contamination.csv"
    per_cell_df.to_csv(csv_path, index=False)
    csv_path = f"{output_dir}/matrix_contamination.csv"
    strength_df.to_csv(csv_path, index=False)
    csv_path = f"{output_dir}/binary_contamination.csv"
    binary_df.to_csv(csv_path, index=False)
    summary["mean_cell_contamination_fraction"] = float(per_cell_df["negative_marker_contamination_fraction"].mean())
    max_contam_val = binary_df.values.max()
    summary["max_type_to_type_contamination_proportion"] = float(max_contam_val)


    #summary
    with open(f"{output_dir}/supervised_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary written to {output_dir}/supervised_summary.json")
    print(f"[INFO] Dataframes generated are stored within respective directories")

    version = subprocess.check_output(
        ["pip", "show", "segtraq"], text=True
    )
    segtraq_version = [l for l in version.splitlines() if l.startswith("Version:")][0].split(": ")[1]

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\n')
        f.write(f'  segtraq: "{segtraq_version}"\n')
        f.write(f'  spatialdata: "{sd.__version__}"\n')
    print("[FINISH] SegTraQ Supervised QC")

if __name__ == "__main__":
    main()
