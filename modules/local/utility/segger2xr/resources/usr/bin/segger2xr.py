#!/usr/bin/env python3
"""
Convert Segger prediction output to XeniumRanger-compatible format.

Reads Segger PREDICT output (transcripts.parquet with segger_cell_id),
produces Baysor-format segmentation CSV, refined transcripts parquet,
and GeoJSON cell boundary polygons for xeniumranger import-segmentation.
"""

import argparse
import json
from pathlib import Path
from typing import List

import pandas as pd
from scipy.spatial import ConvexHull

# Expected columns in transcripts.parquet
REQUIRED_COLUMNS: List[str] = [
    "transcript_id",
    "cell_id",
    "overlaps_nucleus",
    "feature_name",
    "x_location",
    "y_location",
    "z_location",
    "qv",
]

# Column name for segger cell assignment (varies by segger version)
SEGGER_ID_CANDIDATES: List[str] = ["segger_cell_id", "segger_id"]


def refine_transcripts(parquet_path: str) -> pd.DataFrame:
    """
    Read segger PREDICT output and extract cell assignments.
    Supports both 'segger_cell_id' (newer) and 'segger_id' (older) column names.
    """
    parquet_file = Path(parquet_path)
    if not parquet_file.exists():
        raise FileNotFoundError(f"File not found: {parquet_path}")

    df = pd.read_parquet(parquet_file, engine="pyarrow")

    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Find segger cell assignment column
    segger_col = None
    for candidate in SEGGER_ID_CANDIDATES:
        if candidate in df.columns:
            segger_col = candidate
            break
    if segger_col is None:
        raise ValueError(
            f"No segger cell assignment column found. "
            f"Expected one of {SEGGER_ID_CANDIDATES}, got columns: {list(df.columns)}"
        )

    # Replace cell_id with segger assignment
    cell_id_index = df.columns.get_loc("cell_id")
    df = df.drop(columns=["cell_id"])
    segger_series = df.pop(segger_col)
    df.insert(cell_id_index, "cell_id", segger_series)

    return df


def build_cell_map(df: pd.DataFrame, min_transcripts: int = 3) -> dict:
    """
    Build a mapping from raw segger cell IDs to non-numeric string IDs.

    Only includes cells that have:
    - >= min_transcripts assigned transcripts
    - At least one transcript with valid (non-NaN) x/y coordinates

    Cell IDs use "cell-N" format (hyphen + integer) as required by
    xeniumranger's cell ID parser. Non-numeric to avoid polars Int64 inference.
    """
    cell_ids = df["cell_id"].fillna("UNASSIGNED").astype(str)
    is_unassigned = (cell_ids == "UNASSIGNED") | (cell_ids == "") | (cell_ids == "0")
    assigned = cell_ids[~is_unassigned]
    counts = assigned.value_counts()
    enough_tx = set(counts[counts >= min_transcripts].index)

    # Exclude cells with all-NaN coordinates (no spatial info = useless)
    has_coords = df.dropna(subset=["x_location", "y_location"])
    has_coords_ids = set(has_coords["cell_id"].fillna("UNASSIGNED").astype(str))
    valid_cells = sorted(enough_tx & has_coords_ids)

    return {cell: f"cell-{i + 1}" for i, cell in enumerate(valid_cells)}


def to_baysor_csv(df: pd.DataFrame, output_path: str, cell_map: dict) -> None:
    """
    Convert transcript DataFrame to Baysor-compatible CSV format.

    xeniumranger 4.0 import-segmentation --transcript-assignment expects a
    Baysor segmentation CSV with at minimum: transcript_id, cell, is_noise,
    x, y columns. This function maps Xenium/Segger columns to Baysor format.
    """
    baysor_df = pd.DataFrame()
    baysor_df["transcript_id"] = df["transcript_id"]
    baysor_df["x"] = df["x_location"]
    baysor_df["y"] = df["y_location"]
    baysor_df["z"] = df["z_location"]
    baysor_df["gene"] = df["feature_name"]

    cell_ids = df["cell_id"].fillna("UNASSIGNED").astype(str)
    is_unassigned = (cell_ids == "UNASSIGNED") | (cell_ids == "") | (cell_ids == "0")
    baysor_df["cell"] = cell_ids.map(cell_map).fillna("")
    baysor_df["is_noise"] = is_unassigned.astype(int)

    baysor_df.to_csv(output_path, index=False)

    n_assigned = (~is_unassigned).sum()
    n_noise = is_unassigned.sum()
    n_cells = len(cell_map)
    print(
        f"Baysor CSV: {n_assigned} assigned, {n_noise} noise, {n_cells} cells -> {output_path}"
    )


def _make_buffer_polygon(cx: float, cy: float, radius: float = 0.5) -> list:
    """Create a small square polygon around a centroid as fallback."""
    return [
        [cx - radius, cy - radius],
        [cx + radius, cy - radius],
        [cx + radius, cy + radius],
        [cx - radius, cy + radius],
        [cx - radius, cy - radius],  # close ring
    ]


def generate_viz_polygons(df: pd.DataFrame, output_path: str, cell_map: dict) -> None:
    """
    Generate a GeoJSON file with cell boundary polygons.

    Uses ConvexHull when possible; falls back to a small buffer polygon around
    the centroid for cells with < 3 unique points or collinear points.

    Required by xeniumranger import-segmentation when using --transcript-assignment.
    Each feature MUST have a top-level "id" field (xeniumranger reads item["id"]).
    Cell IDs must match those in the Baysor CSV.
    """
    assigned = df[
        df["cell_id"].notna()
        & (df["cell_id"].astype(str) != "UNASSIGNED")
        & (df["cell_id"].astype(str) != "")
    ].copy()

    features = []
    grouped = assigned.groupby("cell_id")

    for cell_id, group in grouped:
        mapped_id = cell_map.get(str(cell_id))
        if mapped_id is None:
            continue

        coords = group[["x_location", "y_location"]].dropna().values

        polygon_coords = None
        if len(coords) >= 3:
            try:
                hull = ConvexHull(coords)
                hull_points = coords[hull.vertices].tolist()
                hull_points.append(hull_points[0])  # close polygon ring
                polygon_coords = hull_points
            except Exception:
                pass

        # Fallback: buffer polygon around centroid
        if polygon_coords is None:
            cx, cy = coords.mean(axis=0).astype(float)
            polygon_coords = _make_buffer_polygon(cx, cy)

        features.append(
            {
                "type": "Feature",
                "id": mapped_id,
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [polygon_coords],
                },
                "properties": {"cell_id": mapped_id},
            }
        )

    geojson = {"type": "FeatureCollection", "features": features}

    with open(output_path, "w") as f:
        json.dump(geojson, f)

    print(f"Generated {len(features)} cell polygons in {output_path}")


def main(input_file: str, prefix: str, min_transcripts: int = 3) -> None:
    """Run the full segger-to-xeniumranger conversion pipeline."""
    Path(prefix).mkdir(parents=True, exist_ok=True)
    transcripts = refine_transcripts(input_file)

    # Build cell ID mapping, filtering cells with < min_transcripts
    cell_map = build_cell_map(transcripts, min_transcripts=min_transcripts)

    # xeniumranger 4.0 expects Baysor-format CSV (not parquet) with is_noise column
    to_baysor_csv(transcripts, f"{prefix}/segmentation.csv", cell_map)

    # Also save the refined parquet for downstream use
    transcripts.to_parquet(f"{prefix}/transcripts.parquet", engine="pyarrow")

    # Generate cell boundary polygons (required companion to --transcript-assignment)
    # Uses ConvexHull when possible; falls back to buffer polygon for edge cases
    generate_viz_polygons(transcripts, f"{prefix}/segmentation_polygons.json", cell_map)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert Segger prediction output to XeniumRanger-compatible format."
    )
    parser.add_argument(
        "--transcripts",
        required=True,
        help="Path to Segger output transcripts parquet file",
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Output directory prefix (sample ID)",
    )
    parser.add_argument(
        "--min-transcripts",
        type=int,
        default=3,
        help="Minimum transcripts per cell (default: 3)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        input_file=args.transcripts,
        prefix=args.prefix,
        min_transcripts=args.min_transcripts,
    )
