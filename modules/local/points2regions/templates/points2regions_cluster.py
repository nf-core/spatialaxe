#!/usr/bin/env python

import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
from points2regions import Points2Regions  # adjust if local or installable


def cluster_points2regions(data, smoothing, num_clusters):
    mdl = Points2Regions(
        data[['X', 'Y']],
        data['gene'],
        pixel_width=1,
        pixel_smoothing=smoothing
    )
    data['clusters'] = mdl.fit_predict(num_clusters=num_clusters, output='marker')
    return data


def plot_clusters(data, smoothing, output_file):
    plt.figure(figsize=(6, 6))
    plt.scatter(
        data['X'],
        data['Y'],
        c=data['clusters'],
        alpha=0.7,
        s=0.5,
        cmap='tab20'
    )
    plt.title(f'Smoothing: {smoothing}')
    plt.axis('off')
    plt.axis('scaled')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)


def main():
    # Inputs passed from Nextflow
    input_csv = "${transcripts}"
    smoothing = int("${smoothing}")
    num_clusters = int("${num_clusters}")

    # Outputs
    clustered_csv = f"clustered_s{smoothing}.csv"
    plot_file = f"cluster_plot_s{smoothing}.png"

    print("[START] Reading transcripts and clustering...")
    data = pd.read_csv(input_csv)
    data = cluster_points2regions(data, smoothing, num_clusters)
    data.to_csv(clustered_csv, index=False)
    print(f"[CLUSTERING DONE] Saved to {clustered_csv}")

    print("[START] Plotting clusters...")
    plot_clusters(data, smoothing, plot_file)
    print(f"[PLOT DONE] Saved to {plot_file}")

    # Version file
    with open("versions.yml", "w") as f:
        f.write(f'"${task.process}":\n')
        f.write('  points2regions: "v1.0.0"\n')

    print("[FINISH]")


if __name__ == "__main__":
    main()
