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
    origin_x: float = 0.0,
    origin_y: float = 0.0,
) -> tuple[float, int]:
    """Return the fraction of transcripts whose (x, y) falls inside a non-zero
    pixel in *mask_tif*, and the count of transcripts dropped for falling
    outside the mask's pixel bounds.

    Reads the mask with tifffile and the transcripts with pandas (lazy chunked
    if the file is very large).

    Transcripts whose pixel coordinates fall outside the mask's bounds are
    excluded from both the numerator and denominator of the fraction, rather
    than clamped to the nearest edge pixel: a transcript off the edge of the
    imaged area carries no information about whether it would have landed in
    a cell, so counting it either way would bias the metric.

    *origin_x* / *origin_y* let the caller specify the transcript-space
    coordinate that corresponds to mask pixel (0, 0), for masks that aren't
    aligned to the transcripts' coordinate origin. They default to (0, 0),
    i.e. "the mask and transcripts already share an origin."
    """
    import pandas as pd
    import tifffile

    mask = tifffile.imread(str(mask_tif))
    h, w = mask.shape[:2]

    inside = 0
    total = 0
    out_of_bounds = 0

    chunksize = 500_000
    for chunk in pd.read_csv(transcripts_csv, usecols=[x_col, y_col], chunksize=chunksize):
        xs = chunk[x_col].values
        ys = chunk[y_col].values

        # Convert micron coordinates to pixel indices, relative to the mask's
        # origin. Cellpose masks are in pixel space; transcripts are in
        # micron space. Divide by pixel_size to convert microns -> pixels.
        # When pixel_size is 1 (default), coords are already in pixels.
        col_f = (xs - origin_x) / pixel_size
        row_f = (ys - origin_y) / pixel_size
        col_idx = np.round(col_f).astype(int)
        row_idx = np.round(row_f).astype(int)

        in_bounds = (col_idx >= 0) & (col_idx < w) & (row_idx >= 0) & (row_idx < h)

        chunk_total = len(chunk)
        chunk_out_of_bounds = int(chunk_total - np.count_nonzero(in_bounds))
        out_of_bounds += chunk_out_of_bounds

        if np.any(in_bounds):
            inside += int(np.count_nonzero(mask[row_idx[in_bounds], col_idx[in_bounds]]))

        # Out-of-bounds transcripts are excluded from the denominator too —
        # we have no information about whether they'd fall in a cell.
        total += int(np.count_nonzero(in_bounds))

    if out_of_bounds:
        print(
            f"WARNING: {out_of_bounds} transcript(s) fell outside the mask's "
            f"pixel bounds and were excluded from the fraction.",
            file=sys.stderr,
        )

    if total == 0:
        return 0.0, out_of_bounds
    return inside / total, out_of_bounds


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
        "--origin-x",
        type=float,
        default=0.0,
        help=(
            "Transcript-space x-coordinate corresponding to mask pixel column 0. "
            "Set this if the mask's pixel origin isn't aligned with the transcripts' "
            "coordinate origin. Default 0.0 (mask and transcripts already aligned)."
        ),
    )
    parser.add_argument(
        "--origin-y",
        type=float,
        default=0.0,
        help=(
            "Transcript-space y-coordinate corresponding to mask pixel row 0. "
            "Set this if the mask's pixel origin isn't aligned with the transcripts' "
            "coordinate origin. Default 0.0 (mask and transcripts already aligned)."
        ),
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

    # --- Validate inputs ---------------------------------------------------
    if args.pixel_size <= 0:
        print(f"ERROR: --pixel-size must be positive, got {args.pixel_size}", file=sys.stderr)
        sys.exit(1)

    # --- Compute fraction -------------------------------------------------
    fraction, out_of_bounds = transcript_in_cell_fraction(
        Path(args.transcripts),
        Path(args.mask),
        x_col=args.x_col,
        y_col=args.y_col,
        pixel_size=args.pixel_size,
        origin_x=args.origin_x,
        origin_y=args.origin_y,
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
            "transcripts_out_of_bounds": out_of_bounds,
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
