#!/usr/bin/env python3
"""
Convert a segmentation mask TIFF to uint32 dtype.

XeniumRanger import-segmentation requires uint32 masks, but upstream
segmenters (e.g. StarDist) often emit int32 labels. This script reads
the input mask, casts it to uint32, and writes the result.
"""

import argparse

import numpy as np
import tifffile


def convert_mask_to_uint32(input_path: str, output_path: str) -> None:
    """
    Read a mask TIFF, cast to uint32, and write to output_path.

    Args:
        input_path: Path to input mask TIFF (any integer dtype).
        output_path: Path where the uint32 mask will be written.
    """
    mask = tifffile.imread(input_path)
    print(f"Input dtype: {mask.dtype}, shape: {mask.shape}, labels: {mask.max()}")
    tifffile.imwrite(output_path, mask.astype(np.uint32))
    print("Output dtype: uint32")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert a segmentation mask TIFF to uint32 dtype."
    )
    parser.add_argument(
        "--input", required=True, help="Path to input mask TIFF"
    )
    parser.add_argument(
        "--output", required=True, help="Path where uint32 mask will be written"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    convert_mask_to_uint32(input_path=args.input, output_path=args.output)
