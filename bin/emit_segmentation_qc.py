#!/usr/bin/env python3
"""Emit segmentation QC JSON from Baysor output.

Reads Baysor segmentation.csv and extracts key metrics:
- Total cells segmented
- Mean/median transcripts per cell
- Transcript assignment rate (assigned / total)

Output: JSON file conforming to the spatialaxe QC schema.
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def compute_segmentation_metrics(
    segmentation_csv: Path,
    transcripts_csv: Path | None = None,
) -> dict:
    """Extract segmentation quality metrics from Baysor output."""
    seg = pd.read_csv(segmentation_csv)

    # Baysor segmentation.csv has a 'cell' column with cell assignments
    # Transcripts with cell=0 are unassigned
    cell_col = "cell"
    if cell_col not in seg.columns:
        # Fallback: try cell_id
        cell_col = "cell_id" if "cell_id" in seg.columns else seg.columns[0]

    total_transcripts = len(seg)
    assigned = seg[cell_col] != 0
    num_assigned = int(assigned.sum())
    assignment_rate = num_assigned / total_transcripts if total_transcripts > 0 else 0.0

    # Per-cell stats (only assigned transcripts)
    cell_counts = seg.loc[assigned, cell_col].value_counts()
    num_cells = len(cell_counts)
    mean_transcripts_per_cell = float(cell_counts.mean()) if num_cells > 0 else 0.0
    median_transcripts_per_cell = float(cell_counts.median()) if num_cells > 0 else 0.0

    return {
        "total_transcripts": total_transcripts,
        "num_cells": num_cells,
        "assigned_transcripts": num_assigned,
        "assignment_rate": round(assignment_rate, 4),
        "mean_transcripts_per_cell": round(mean_transcripts_per_cell, 1),
        "median_transcripts_per_cell": round(median_transcripts_per_cell, 1),
    }


def main():
    parser = argparse.ArgumentParser(description="Emit segmentation QC JSON.")
    parser.add_argument("--segmentation", required=True, help="Baysor segmentation.csv")
    parser.add_argument("--sample-id", required=True, help="Sample identifier")
    parser.add_argument("--output", default="segmentation_qc.json", help="Output JSON path")
    args = parser.parse_args()

    metrics = compute_segmentation_metrics(Path(args.segmentation))

    qc_record = {
        "sample_id": args.sample_id,
        "stage": "segmentation",
        "metrics": metrics,
    }

    Path(args.output).write_text(json.dumps(qc_record, indent=2) + "\n")
    print(f"Segmentation QC: {metrics['num_cells']} cells, "
          f"assignment rate {metrics['assignment_rate']:.2%}")


if __name__ == "__main__":
    main()
