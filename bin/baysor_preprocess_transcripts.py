#!/usr/bin/env python3
"""
Preprocess Xenium transcripts for Baysor segmentation.

Filters transcripts based on quality score and spatial coordinate thresholds,
removes negative control probes, and outputs filtered CSV for Baysor compatibility.
"""

import argparse
import os

import pandas as pd


def filter_transcripts(
    transcripts: str,
    min_qv: float = 20.0,
    min_x: float = 0.0,
    max_x: float = 24000.0,
    min_y: float = 0.0,
    max_y: float = 24000.0,
    prefix: str = "",
) -> None:
    """
    Filter transcripts based on the specified thresholds.

    Args:
        transcripts: Path to transcripts parquet file
        min_qv: Minimum Q-Score to pass filtering
        min_x: Minimum x-coordinate threshold
        max_x: Maximum x-coordinate threshold
        min_y: Minimum y-coordinate threshold
        max_y: Maximum y-coordinate threshold
        prefix: Output directory prefix
    """
    df = pd.read_parquet(transcripts, engine="pyarrow")

    # filter transcripts df with thresholds, ignore negative controls
    filtered_df = df[
        (df["qv"] >= min_qv)
        & (df["x_location"] >= min_x)
        & (df["x_location"] <= max_x)
        & (df["y_location"] >= min_y)
        & (df["y_location"] <= max_y)
        & (~df["feature_name"].str.startswith("NegControlProbe_"))
        & (~df["feature_name"].str.startswith("antisense_"))
        & (~df["feature_name"].str.startswith("NegControlCodeword_"))
        & (~df["feature_name"].str.startswith("BLANK_"))
    ]

    # change cell_id of cell-free transcripts to "0" (Baysor's no-cell sentinel).
    # Modern Xenium stores cell_id as a string ("UNASSIGNED" for cell-free transcripts);
    # legacy Xenium used integer -1. Normalize to string and handle both cases — pandas 3
    # rejects mixing int values into a string-dtype column.
    filtered_df["cell_id"] = filtered_df["cell_id"].astype(str)
    neg_cell_row = filtered_df["cell_id"].isin(["-1", "UNASSIGNED"])
    filtered_df.loc[neg_cell_row, "cell_id"] = "0"

    # Output filtered transcripts as CSV for Baysor 0.7.1 compatibility.
    # Baysor's Julia Parquet.jl cannot read modern pyarrow Parquet files
    # (pyarrow 15+ writes size_statistics Thrift field 16 unconditionally,
    # which Baysor's old Thrift deserializer doesn't recognize).
    os.makedirs(prefix, exist_ok=True)
    filtered_df.to_csv(f"{prefix}/filtered_transcripts.csv", index=False)

    return None


def main() -> None:
    """
    Run preprocess transcripts as nf module.
    """
    parser = argparse.ArgumentParser(
        description="Preprocess Xenium transcripts for Baysor"
    )
    parser.add_argument(
        "--transcripts", required=True, help="Path to transcripts parquet file"
    )
    parser.add_argument("--prefix", required=True, help="Output directory prefix")
    parser.add_argument(
        "--min-qv",
        type=float,
        default=20.0,
        help="Minimum Q-Score threshold (default: 20.0)",
    )
    parser.add_argument(
        "--min-x",
        type=float,
        default=0.0,
        help="Minimum x-coordinate threshold (default: 0.0)",
    )
    parser.add_argument(
        "--max-x",
        type=float,
        default=24000.0,
        help="Maximum x-coordinate threshold (default: 24000.0)",
    )
    parser.add_argument(
        "--min-y",
        type=float,
        default=0.0,
        help="Minimum y-coordinate threshold (default: 0.0)",
    )
    parser.add_argument(
        "--max-y",
        type=float,
        default=24000.0,
        help="Maximum y-coordinate threshold (default: 24000.0)",
    )
    args = parser.parse_args()

    filter_transcripts(
        transcripts=args.transcripts,
        min_qv=args.min_qv,
        min_x=args.min_x,
        max_x=args.max_x,
        min_y=args.min_y,
        max_y=args.max_y,
        prefix=args.prefix,
    )

    return None


if __name__ == "__main__":
    main()
