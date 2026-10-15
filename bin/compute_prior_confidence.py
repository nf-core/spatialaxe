#!/usr/bin/env python3
"""Compute prior segmentation confidence for Baysor from transcript-in-cell fraction.

Decision 1C: Uses transcript-in-cell fraction — the proportion of transcripts
that fall inside a Cellpose segmentation mask — as the quality metric.

Decision 2C: Mapping function is user-configurable via --mapping; ships with
a sensible linear default.  Users can switch to 'binned' or supply a custom
two-column CSV (metric_value,confidence).

Decision 3A: The auto-computed value is the default.  If the user hard-sets
`baysor_prior_confidence` in the pipeline params, this script is bypassed.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def transcript_in_cell_fraction(
    transcripts_csv: Path,
    mask_tif: Path,
    x_col: str = "x_location",
    y_col: str = "y_location",
    pixel_size: float = 1.0,
) -> float:
    """Return the fraction of transcripts whose (x, y) falls inside a non-zero
    pixel in *mask_tif*.

    Reads the mask with tifffile and the transcripts with pandas (lazy chunked
    if the file is very large).
    """
    import pandas as pd
    import tifffile

    mask = tifffile.imread(str(mask_tif))
    h, w = mask.shape[:2]

    inside = 0
    total = 0

    chunksize = 500_000
    for chunk in pd.read_csv(transcripts_csv, usecols=[x_col, y_col], chunksize=chunksize):
        xs = chunk[x_col].values
        ys = chunk[y_col].values

        # Convert micron coordinates to pixel indices.
        # Cellpose masks are in pixel space; transcripts are in micron space.
        # Divide by pixel_size to convert microns → pixels.
        # When pixel_size is 1 (default), coords are already in pixels.
        col_idx = np.clip(np.round(xs / pixel_size).astype(int), 0, w - 1)
        row_idx = np.clip(np.round(ys / pixel_size).astype(int), 0, h - 1)

        inside += int(np.count_nonzero(mask[row_idx, col_idx]))
        total += len(chunk)

    if total == 0:
        return 0.0
    return inside / total


def map_linear(fraction: float, floor: float = 0.05, ceiling: float = 0.95) -> float:
    """Rescale fraction ∈ [0, 1] to confidence ∈ [floor, ceiling]."""
    return float(np.clip(floor + fraction * (ceiling - floor), floor, ceiling))


def map_binned(fraction: float) -> float:
    """Map fraction to low / medium / high bins."""
    if fraction < 0.3:
        return 0.2
    elif fraction < 0.7:
        return 0.5
    else:
        return 0.8


def map_custom(fraction: float, csv_path: Path) -> float:
    """Linearly interpolate from a user-supplied two-column CSV
    (metric_value,confidence).
    """
    import pandas as pd

    lut = pd.read_csv(csv_path).sort_values("metric_value")
    return float(np.interp(fraction, lut["metric_value"].values, lut["confidence"].values))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute Baysor prior_segmentation_confidence from transcript-in-cell fraction."
    )
    parser.add_argument("--transcripts", required=True, help="Path to preprocessed transcripts CSV.")
    parser.add_argument("--mask", required=True, help="Path to Cellpose segmentation mask (TIF).")
    parser.add_argument(
        "--mapping",
        default="linear",
        help=(
            "Mapping strategy: 'linear' (default), 'binned', or path to a "
            "two-column CSV (metric_value,confidence) for custom interpolation."
        ),
    )
    parser.add_argument(
        "--pixel-size",
        type=float,
        default=1.0,
        help="Microns per pixel.  Transcript coords are divided by this before lookup.",
    )
    parser.add_argument(
        "--x-col", default="x_location", help="Column name for x coordinates."
    )
    parser.add_argument(
        "--y-col", default="y_location", help="Column name for y coordinates."
    )
    parser.add_argument(
        "--sample-id", default="unknown", help="Sample identifier for QC reporting."
    )
    parser.add_argument(
        "--output", default="prior_confidence.json", help="Output JSON path."
    )
    args = parser.parse_args()

    # --- Compute fraction -------------------------------------------------
    fraction = transcript_in_cell_fraction(
        Path(args.transcripts),
        Path(args.mask),
        x_col=args.x_col,
        y_col=args.y_col,
        pixel_size=args.pixel_size,
    )

    # --- Map to confidence ------------------------------------------------
    mapping = args.mapping
    if mapping == "linear":
        confidence = map_linear(fraction)
    elif mapping == "binned":
        confidence = map_binned(fraction)
    elif Path(mapping).is_file():
        confidence = map_custom(fraction, Path(mapping))
    else:
        print(f"ERROR: unknown mapping '{mapping}'; expected 'linear', 'binned', or a CSV path.", file=sys.stderr)
        sys.exit(1)

    # --- Write output -----------------------------------------------------
    result = {
        "sample_id": args.sample_id,
        "stage": "prior_confidence",
        "metrics": {
            "transcript_in_cell_fraction": round(fraction, 4),
            "prior_segmentation_confidence": round(confidence, 4),
        },
        "mapping": mapping if mapping in ("linear", "binned") else "custom",
    }
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n")

    print(
        f"Transcript-in-cell fraction: {fraction:.4f} → "
        f"prior_segmentation_confidence: {confidence:.4f} (mapping={result['mapping']})"
    )


if __name__ == "__main__":
    main()
