#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

def plot_clusters(data, smoothing, output_file="cluster_plot.png"):
    """
    Plot spatial clusters from Points2Regions output.

    Args:
        data (pd.DataFrame): Data with 'X', 'Y', and 'Clusters' columns.
        smoothing (int): Smoothing value used (for the plot title).
        output_file (str): Path to save the plot image.
    """
    plt.figure(figsize=(6, 6))
    plt.scatter(
        data['X'],
        data['Y'],
        c=data['Clusters'],
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
    input_csv = "${clustered}"
    smoothing = int("${smoothing}")
    output_file = "cluster_plot.png"

    data = pd.read_csv(input_csv)
    plot_clusters(data, smoothing, output_file)

    # Save version info
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('  points2regions_plot: "v1.0.0"\\n')

if __name__ == "__main__":
    main()
