#!/usr/bin/env python3
"""Merge two spatialdata bundles to create a layered spatialdata object."""

import argparse
import json
import os
import shutil

import spatialdata


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Merge two spatialdata bundles")
    parser.add_argument("--raw-bundle", required=True, help="Path to raw spatialdata bundle")
    parser.add_argument("--redefined-bundle", required=True, help="Path to redefined spatialdata bundle")
    parser.add_argument("--prefix", required=True, help="Output prefix (sample ID)")
    parser.add_argument("--output-folder", required=True, help="Output folder name")
    return parser.parse_args()


def main():
    """Run spatialdata merge."""
    args = parse_args()
    print("[START]")

    output_dir = f"spatialdata/{args.prefix}/{args.output_folder}"

    # Ensure the output folder exists
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    # Copy the entire reference bundle as is
    for root, _, files in os.walk(args.raw_bundle):
        rel_path = os.path.relpath(root, args.raw_bundle)
        target_path = os.path.join(output_dir, rel_path)
        os.makedirs(target_path, exist_ok=True)
        for file in files:
            shutil.copy(os.path.join(root, file), os.path.join(target_path, file))

    # Rename folders in Points, Shapes, and Tables to raw_*
    for category in ["points", "shapes", "tables"]:
        category_path = os.path.join(output_dir, category)
        if os.path.exists(category_path):
            for folder in next(os.walk(category_path))[1]:
                old_path = os.path.join(category_path, folder)
                print(folder)
                new_path = os.path.join(category_path, f"raw_{folder}")
                os.rename(old_path, new_path)

    # Copy folders from redefined_bundle and rename them as redefined_*
    for category in ["points", "shapes", "tables"]:
        add_category_path = os.path.join(args.redefined_bundle, category)
        output_category_path = os.path.join(output_dir, category)
        os.makedirs(output_category_path, exist_ok=True)

        if os.path.exists(add_category_path):
            for folder in next(os.walk(add_category_path))[1]:
                src_folder = os.path.join(add_category_path, folder)
                dest_folder = os.path.join(output_category_path, f"redefined_{folder}")
                shutil.copytree(src_folder, dest_folder)

    # Invalidate consolidated metadata in zarr.json -- the directory renames above
    # made the element paths in the metadata stale (e.g., 'points/transcripts' ->
    # 'points/raw_transcripts'). Without consolidated metadata, sd.read_zarr()
    # discovers elements by scanning the filesystem directly.
    zarr_json = os.path.join(output_dir, "zarr.json")
    if os.path.exists(zarr_json):
        with open(zarr_json) as f:
            meta = json.load(f)
        if "consolidated_metadata" in meta:
            del meta["consolidated_metadata"]
            with open(zarr_json, "w") as f:
                json.dump(meta, f)
            print("[NOTE] Removed stale consolidated metadata from zarr.json")

    print("[FINISH]")


if __name__ == "__main__":
    main()
