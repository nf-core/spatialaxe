#!/usr/bin/env python

import sys
import os
import pandas as pd
from points2regions import Points2Regions  # adjust if function is local

def cluster_points2regions():
    print("[START]")

    input_csv = "${transcripts}"
    output_csv = "clustered.csv"
    smoothing = int("${smoothing}")
    num_clusters = int("${num_clusters}")

    # Read input
    data = pd.read_csv(input_csv)

    # Run clustering
    mdl = Points2Regions(
        data[['X', 'Y']],
        data['gene'],
        pixel_width=1,
        pixel_smoothing=smoothing
    )

    data['clusters'] = mdl.fit_predict(num_clusters=num_clusters, output='marker')

    # Save result
    data.to_csv(output_csv, index=False)

    # Write version info
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('  points2regions_cluster: "v1.0.0"\\n')

    print("[FINISH]")

if __name__ == "__main__":
    cluster_points2regions()
