#!/usr/bin/env python3
"""
Resize a segmentation TIFF mask to match transcript coordinates.

This script rescales a segmentation mask image to match the coordinate
space of Xenium transcript data using microns-per-pixel metadata.
"""

import argparse
import json
import os
from typing import Tuple

import numpy as np
import pandas as pd
import tifffile
from skimage.transform import resize


def read_mask(mask_path: str) -> np.ndarray:
    """Read the segmentation mask from a TIFF file."""
    print(f"Reading mask: {mask_path}")
    mask = tifffile.imread(mask_path)
    print(f"Mask shape: {mask.shape}, dtype: {mask.dtype}")
    return mask


def read_transcript_bounds(transcript_path: str) -> Tuple[float, float, float, float]:
    """Read transcript coordinates and return their bounding box."""
    print(f"Reading transcripts: {transcript_path}")
    if transcript_path.endswith(".parquet"):
        transcripts = pd.read_parquet(transcript_path, columns=["x_location", "y_location"])
    else:
        transcripts = pd.read_csv(transcript_path)

    if "x_location" not in transcripts.columns or "y_location" not in transcripts.columns:
        raise ValueError("Transcript file must contain 'x_location' and 'y_location' columns.")

    x_min, x_max = transcripts["x_location"].min(), transcripts["x_location"].max()
    y_min, y_max = transcripts["y_location"].min(), transcripts["y_location"].max()

    print(f"Transcript bounds: X=({x_min:.2f}, {x_max:.2f}), Y=({y_min:.2f}, {y_max:.2f})")
    return x_min, x_max, y_min, y_max


def read_microns_per_pixel(metadata_path: str) -> float:
    """Extract microns_per_pixel or pixel_size from metadata JSON."""
    print(f"Reading metadata: {metadata_path}")
    with open(metadata_path, "r") as f:
        metadata = json.load(f)

    mpp = metadata.get("microns_per_pixel") or metadata.get("pixel_size")
    if mpp is None:
        raise KeyError("Metadata JSON must contain 'microns_per_pixel' or 'pixel_size'.")

    print(f"Microns per pixel: {mpp}")
    return float(mpp)


def compute_target_size(
    x_min: float, x_max: float, y_min: float, y_max: float, microns_per_pixel: float
) -> Tuple[int, int]:
    """Compute new image size (in pixels) to cover given coordinates."""
    new_width = int(round((x_max - x_min) / microns_per_pixel))
    new_height = int(round((y_max - y_min) / microns_per_pixel))
    print(f"Target image size: {new_width} x {new_height} pixels")
    return new_height, new_width


def resize_mask(mask: np.ndarray, new_shape: Tuple[int, int]) -> np.ndarray:
    """Resize mask using nearest-neighbor interpolation (preserve labels)."""
    print("Resizing mask...")
    resized = resize(
        mask,
        new_shape,
        order=0,  # nearest neighbor to preserve segmentation labels
        preserve_range=True,
        anti_aliasing=False,
    ).astype(mask.dtype)
    print(f"Resized shape: {resized.shape}")
    return resized


def main(mask_path: str, transcripts_path: str, metadata_path: str, output_path: str) -> None:
    """Resize segmentation mask to match Xenium coordinate space."""
    # Validate input files
    for path in [mask_path, transcripts_path, metadata_path]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")

    # Load data
    mask = read_mask(mask_path)
    x_min, x_max, y_min, y_max = read_transcript_bounds(transcripts_path)
    microns_per_pixel = read_microns_per_pixel(metadata_path)

    # Compute physical mask size
    height, width = mask.shape
    print(f"Original mask size: {width * microns_per_pixel:.2f} x {height * microns_per_pixel:.2f} um")

    # Compute target size
    new_height, new_width = compute_target_size(x_min, x_max, y_min, y_max, microns_per_pixel)

    # Resize and save
    resized_mask = resize_mask(mask, (new_height, new_width))
    tifffile.imwrite(output_path, resized_mask)

    print(f"Saved resized mask -> {output_path}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Resize a segmentation TIFF mask to match transcript coordinates."
    )
    parser.add_argument("--mask", required=True, help="Path to segmentation mask TIFF")
    parser.add_argument("--transcripts", required=True, help="Path to transcripts file")
    parser.add_argument("--metadata", required=True, help="Path to metadata JSON")
    parser.add_argument("--prefix", required=True, help="Output directory prefix")
    parser.add_argument("--mask-filename", required=True, help="Original mask filename for output naming")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    os.makedirs(args.prefix, exist_ok=True)
    output_mask: str = os.path.join(args.prefix, f"resized_{args.mask_filename}.tif")

    main(
        mask_path=args.mask,
        transcripts_path=args.transcripts,
        metadata_path=args.metadata,
        output_path=output_mask,
    )
