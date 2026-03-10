#!/usr/bin/env python

"""Merge two spatialdata bundles to create a layered spatialdata object."""

import json
import os
import shutil
import spatialdata

def main():

    print("[START]")

    raw_bundle = "${raw_bundle}"
    redefined_bundle = "${redefined_bundle}"
    output_path = "${prefix}"
    output_folder = "${outputfolder}"

    # Ensure the output folder exists
    if os.path.exists(f"spatialdata/{output_path}/{output_folder}"):
        shutil.rmtree(f"spatialdata/{output_path}/{output_folder}")
    os.makedirs(f"spatialdata/{output_path}/{output_folder}")

    # Copy the entire reference bundle as is
    for root, _, files in os.walk(raw_bundle):
        rel_path = os.path.relpath(root, raw_bundle)
        target_path = os.path.join(f"spatialdata/{output_path}/{output_folder}", rel_path)
        os.makedirs(target_path, exist_ok=True)
        for file in files:
            shutil.copy(os.path.join(root, file), os.path.join(target_path, file))

    # Rename folders in Points, Shapes, and Tables to raw_*
    for category in ["points", "shapes", "tables"]:
        category_path = os.path.join(f"spatialdata/{output_path}/{output_folder}", category)
        if os.path.exists(category_path):
            for folder in next(os.walk(category_path))[1]: #os.listdir(category_path):
                old_path = os.path.join(category_path, folder)
                print(folder)
                new_path = os.path.join(category_path, f"raw_{folder}")
                os.rename(old_path, new_path)

    # Copy folders from redefined_bundle and rename them as redefined_*
    for category in ["points", "shapes", "tables"]:
        add_category_path = os.path.join(redefined_bundle, category)
        output_category_path = os.path.join(f"spatialdata/{output_path}/{output_folder}", category)
        os.makedirs(output_category_path, exist_ok=True)

        if os.path.exists(add_category_path):
            for folder in next(os.walk(add_category_path))[1]:
                src_folder = os.path.join(add_category_path, folder)
                dest_folder = os.path.join(output_category_path, f"redefined_{folder}")
                shutil.copytree(src_folder, dest_folder)

    # Invalidate consolidated metadata in zarr.json — the directory renames above
    # made the element paths in the metadata stale (e.g., 'points/transcripts' →
    # 'points/raw_transcripts'). Without consolidated metadata, sd.read_zarr()
    # discovers elements by scanning the filesystem directly.
    output_dir = f"spatialdata/{output_path}/{output_folder}"
    zarr_json = os.path.join(output_dir, "zarr.json")
    if os.path.exists(zarr_json):
        with open(zarr_json) as f:
            meta = json.load(f)
        if "consolidated_metadata" in meta:
            del meta["consolidated_metadata"]
            with open(zarr_json, "w") as f:
                json.dump(meta, f)
            print("[NOTE] Removed stale consolidated metadata from zarr.json")

    #Output version information
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'spatialdata: "{spatialdata.__version__}"\\n')

    print("[FINISH]")

if __name__ == "__main__":
    main()
