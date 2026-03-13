#!/usr/bin/env python

"""Compute baseline statistics on spatialdata object for QC."""

import os
import segtraq
import spatialdata as sd
import json
import subprocess

def main():
    print("[START] SegTraQ Baseline QC")
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

    print(f"[INFO] Computing baseline QC metrics")
    summary = {}

    #number of cells
    n_cells = st.bl.num_cells()
    summary["num_cells"] = int(n_cells)
    print(f"  num_cells: {n_cells}")

    #number of transcripts
    n_transcripts = st.bl.num_transcripts()
    summary["num_transcripts"] = int(n_transcripts)
    print(f" num_transcripts: {n_transcripts}")

    #number of genes
    n_genes = st.bl.num_genes()
    summary["num_genes"] = int(n_genes)
    print(f" num_genes: {n_genes}")

    #percentage of assigned transcripts
    percentage_unassgn_transcripts = st.bl.perc_unassigned_transcripts()
    summary["percent_unassigned_transcripts"] = int(percentage_unassgn_transcripts)
    print(f" percent_unassigned_transcripts: {percentage_unassgn_transcripts}")

    #unassigned transcripts per gene
    unassgn_transcripts_per_gene = st.bl.perc_unassigned_transcripts_per_gene()
    summary["unassigned_transcripts_per_gene"] = int(unassgn_transcripts_per_gene)
    print(f" unassigned_transcripts_per_gene: {unassgn_transcripts_per_gene}")

    #transcripts per cell
    transcripts_per_cell = st.bl.transcripts_per_cell()
    summary["transcripts_per_cell"] = int(transcripts_per_cell)
    print(f" transcripts_per_cell: {transcripts_per_cell}")

    #genes per celll
    genes_per_cell = st.bl.genes_per_cell()
    summary["genes_per_cell"] = int(genes_per_cell)
    print(f" genes_per_cell: {genes_per_cell}")

    #transcript density
    transcript_density = st.bl.transcript_density()
    summary["transcript_density"] = int(transcript_density)
    print(f" transcript_density: {transcript_density}")

    #mean transcripts per gene cell
    mean_transcripts_per_gene_per_cell = st.bl.mean_transcripts_per_gene_per_cell()
    summary["mean_transcripts_per_gene_per_cell"] = int(mean_transcripts_per_gene_per_cell)
    print(f" mean_transcripts_per_gene_per_cell: {mean_transcripts_per_gene_per_cell}")

    #morphological features
    morpho_features = st.bl.morphological_features()
    summary["morpho_features"] = int(morpho_features)
    print(f" morpho_features: {morpho_features}")

    #summary
    with open(f"{output_dir}/baseline_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary written to {output_dir}/baseline_summary.json")

    version = subprocess.check_output(
        ["pip", "show", "segtraq"], text=True
    )
    segtraq_version = [l for l in version.split("\\n") if l.startswith("Version:")][0].split(": ")[1]

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'  segtraq: "{segtraq_version}"\\n')
        f.write(f'  spatialdata: "{sd.__version__}"\\n')
    print("[FINISH] SegTraQ Baseline QC")

if __name__ == "__main__":
    main()










