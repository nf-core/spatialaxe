#!/usr/bin/env python3
"""
Extract a single channel (e.g., DAPI) from a multi-channel OME-TIFF.

Xenium morphology_focus.ome.tif has multiple channels (DAPI, boundary,
interior). Single-channel segmenters such as StarDist 2D_versatile_fluo
expect one channel as input. This script reads the input image, slices
the requested channel, and writes the result.
"""

import argparse

import tifffile


def extract_channel(input_path: str, output_path: str, channel_index: int) -> None:
    """
    Read an OME-TIFF, extract a single channel, and write the result.

    Args:
        input_path: Path to multi-channel OME-TIFF morphology image.
        output_path: Path where the single-channel TIFF will be written.
        channel_index: Index of the channel to extract.
    """
    img = tifffile.imread(input_path)
    orig_shape = img.shape

    if img.ndim == 3:
        img = img[channel_index]
    elif img.ndim == 4:
        img = img[0, channel_index]

    tifffile.imwrite(output_path, img)
    print(f"Input shape: {orig_shape} -> extracted channel {channel_index}: {img.shape}")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract a single channel from a multi-channel OME-TIFF."
    )
    parser.add_argument(
        "--input", required=True, help="Path to multi-channel OME-TIFF morphology image"
    )
    parser.add_argument(
        "--output", required=True, help="Path where the single-channel TIFF will be written"
    )
    parser.add_argument(
        "--channel-index", type=int, default=0, help="Channel index to extract (default: 0)"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    extract_channel(
        input_path=args.input,
        output_path=args.output,
        channel_index=args.channel_index,
    )
