#!/usr/bin/env python3
"""Emit Xenium bundle read QC JSON.

Reads metrics_summary.csv and experiment.xenium from a Xenium bundle to
extract key read-level and experiment-level metrics.

Output: JSON file conforming to the spatialaxe QC schema.
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def extract_bundle_metrics(bundle_path: Path) -> dict:
    """Extract QC metrics from a Xenium bundle."""
    metrics = {}

    # metrics_summary.csv
    metrics_csv = bundle_path / "metrics_summary.csv"
    if metrics_csv.exists():
        df = pd.read_csv(metrics_csv)
        # metrics_summary.csv typically has columns: Category, Metric, Value
        # or is a wide-format single-row table
        if "Metric" in df.columns and "Value" in df.columns:
            for _, row in df.iterrows():
                key = str(row["Metric"]).strip().lower().replace(" ", "_").replace("-", "_")
                try:
                    metrics[key] = float(row["Value"])
                except (ValueError, TypeError):
                    metrics[key] = str(row["Value"])
        else:
            # Wide format: column names are metric names
            for col in df.columns:
                key = col.strip().lower().replace(" ", "_").replace("-", "_")
                val = df[col].iloc[0]
                try:
                    metrics[key] = float(val)
                except (ValueError, TypeError):
                    metrics[key] = str(val)

    # experiment.xenium for metadata
    exp_file = bundle_path / "experiment.xenium"
    if exp_file.exists():
        try:
            exp = json.loads(exp_file.read_text())
            if "num_cells" in exp:
                metrics["xenium_num_cells"] = int(exp["num_cells"])
            if "panel_tissue_type" in exp:
                metrics["tissue_type"] = str(exp["panel_tissue_type"])
        except (json.JSONDecodeError, OSError):
            pass

    # transcripts.parquet row count (if available, without loading all data)
    transcripts_pq = bundle_path / "transcripts.parquet"
    if transcripts_pq.exists():
        try:
            import pyarrow.parquet as pq
            pf = pq.ParquetFile(str(transcripts_pq))
            metrics["total_transcripts"] = pf.metadata.num_rows
        except (ImportError, OSError):
            pass

    return metrics


def main():
    parser = argparse.ArgumentParser(description="Emit Xenium bundle read QC JSON.")
    parser.add_argument("--bundle", required=True, help="Path to Xenium bundle directory")
    parser.add_argument("--sample-id", required=True, help="Sample identifier")
    parser.add_argument("--stage", default="bundle_qc", help="Stage name (default: bundle_qc)")
    parser.add_argument("--output", default="bundle_qc.json", help="Output JSON path")
    args = parser.parse_args()

    metrics = extract_bundle_metrics(Path(args.bundle))

    qc_record = {
        "sample_id": args.sample_id,
        "stage": args.stage,
        "metrics": metrics,
    }

    Path(args.output).write_text(json.dumps(qc_record, indent=2) + "\n")
    n_metrics = len(metrics)
    print(f"Bundle QC: extracted {n_metrics} metrics from {args.bundle}")


if __name__ == "__main__":
    main()
