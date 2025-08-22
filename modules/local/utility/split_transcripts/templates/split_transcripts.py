#!/usr/bin/env python3

import pandas as pd
from typing import List


def compute_quantile_ranges(df: pd.DataFrame, col: str, n_bins: int) -> List:
    """
    Compute the bin edges for `df[col]` such that each of the n_bins
    has ~equal count of points. Returns a list of (min, max) tuples.
    """
    _, bins = pd.qcut(df[col], q=n_bins, retbins=True, duplicates='drop')

    ranges = [(bins[i], bins[i+1]) for i in range(len(bins)-1)]

    return ranges


def make_tiles(df: pd.DataFrame, x_bins: int, y_bins: int) -> pd.DataFrame:
    """
    Produce a DataFrame with one row per tile:
      tile_id, x_min, x_max, y_min, y_max
    """
    x_ranges = compute_quantile_ranges(df, 'x_location', x_bins)
    y_ranges = compute_quantile_ranges(df, 'y_location', y_bins)

    tiles = []
    for ix, (x_min, x_max) in enumerate(x_ranges, start=1):
        for iy, (y_min, y_max) in enumerate(y_ranges, start=1):
            tiles.append({
                'tile_id': f'{ix}_{iy}',
                'x_min': x_min,
                'x_max': x_max,
                'y_min': y_min,
                'y_max': y_max
            })

    return pd.DataFrame(tiles)


def generate_version_yml() -> None:
    with open("versions.yml", "w") as yml:
        yml.write('"${task.process}":\\n')
        yml.write("Baysor-Split Transcripts: 1.0.0'\\n")

    return None


def main(
    transcripts: str,
    x_bins: int = 10,
    y_bins: int = 10
) -> None:
    """
    Generate splits
    """
    # read parquet file
    df = pd.read_parquet(transcripts, engine='fastparquet')

    # compute tiles
    tiles_df = make_tiles(df, x_bins, y_bins)

    # save parquet file
    tiles_df.to_csv("splits.csv", index=False)

    # generate version yml
    generate_version_yml()

    return None


if __name__ == "__main__":

    transcripts: str = "${transcripts}"
    x_bins: int = "${x_bins}"
    y_bins: int = "${y_bins}"

    main(transcripts=transcripts, x_bins=x_bins, y_bins=y_bins)
