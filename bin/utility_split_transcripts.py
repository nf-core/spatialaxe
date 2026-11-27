#!/usr/bin/env python3
"""
Split transcript coordinates into spatial tiles.

Reads a Xenium transcripts.parquet file and computes quantile-based spatial
tiles, writing a splits.csv with tile boundaries.
"""

import argparse
import os
from typing import List

import pandas as pd


def compute_quantile_ranges(df: pd.DataFrame, col: str, n_bins: int) -> List:
    """
    Compute the bin edges for `df[col]` such that each of the n_bins
    has ~equal count of points. Returns a list of (min, max) tuples.
    """
    _, bins = pd.qcut(df[col], q=n_bins, retbins=True, duplicates="drop")

    ranges = [(bins[i], bins[i + 1]) for i in range(len(bins) - 1)]

    return ranges


def make_tiles(df: pd.DataFrame, x_bins: int, y_bins: int) -> pd.DataFrame:
    """
    Produce a DataFrame with one row per tile:
      tile_id, x_min, x_max, y_min, y_max
    """
    x_ranges = compute_quantile_ranges(df, "x_location", x_bins)
    y_ranges = compute_quantile_ranges(df, "y_location", y_bins)

    tiles = []
    for ix, (x_min, x_max) in enumerate(x_ranges, start=1):
        for iy, (y_min, y_max) in enumerate(y_ranges, start=1):
            tiles.append(
                {
                    "tile_id": f"{ix}_{iy}",
                    "x_min": x_min,
                    "x_max": x_max,
                    "y_min": y_min,
                    "y_max": y_max,
                }
            )

    return pd.DataFrame(tiles)


def main(
    transcripts: str,
    x_bins: int = 10,
    y_bins: int = 10,
    prefix: str = "",
) -> None:
    """Generate spatial tile splits from transcript coordinates."""
    # read parquet file
    df = pd.read_parquet(transcripts, engine="fastparquet")

    # compute tiles
    tiles_df = make_tiles(df, x_bins, y_bins)

    # save csv file
    os.makedirs(prefix, exist_ok=True)
    tiles_df.to_csv(f"{prefix}/splits.csv", index=False)

    return None


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Split transcript coordinates into spatial tiles."
    )
    parser.add_argument(
        "--transcripts",
        required=True,
        help="Path to transcripts parquet file",
    )
    parser.add_argument(
        "--x-bins",
        type=int,
        required=True,
        help="Number of bins along X axis",
    )
    parser.add_argument(
        "--y-bins",
        type=int,
        required=True,
        help="Number of bins along Y axis",
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Output directory prefix",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        transcripts=args.transcripts,
        x_bins=args.x_bins,
        y_bins=args.y_bins,
        prefix=args.prefix,
    )
