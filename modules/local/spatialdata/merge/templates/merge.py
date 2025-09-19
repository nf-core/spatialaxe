#!/usr/bin/env python

"""Merge two spatialdata bundles to create a layered spatialdata object."""

import os
import shutil
import spatialdata

def main():

    print("[START]")

    reference_bundle = "${ref_bundle}"
    add_bundle = "${add_bundle}"
    output_folder = "./spatialdata_merged"

    # Ensure the output folder exists
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    os.makedirs(output_folder)

    # Copy the entire reference bundle as is
    for root, _, files in os.walk(reference_bundle):
        rel_path = os.path.relpath(root, reference_bundle)
        target_path = os.path.join(output_folder, rel_path)
        os.makedirs(target_path, exist_ok=True)
        for file in files:
            shutil.copy(os.path.join(root, file), os.path.join(target_path, file))

    # Rename folders in Points, Shapes, and Tables to raw_*
    for category in ["points", "shapes", "tables"]:
        category_path = os.path.join(output_folder, category)
        if os.path.exists(category_path):
            for folder in next(os.walk(category_path))[1]: #os.listdir(category_path):
                old_path = os.path.join(category_path, folder)
                print(folder)
                new_path = os.path.join(category_path, f"raw_{folder}")
                os.rename(old_path, new_path)

    # Copy folders from add_bundle and rename them as reference_*
    for category in ["points", "shapes", "tables"]:
        add_category_path = os.path.join(add_bundle, category)
        output_category_path = os.path.join(output_folder, category)
        os.makedirs(output_category_path, exist_ok=True)

        if os.path.exists(add_category_path):
            for folder in next(os.walk(add_category_path))[1]:
                src_folder = os.path.join(add_category_path, folder)
                dest_folder = os.path.join(output_category_path, f"reference_{folder}")
                shutil.copytree(src_folder, dest_folder)

    #Output version information
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'spatialdata: "{spatialdata.__version__}"\\n')

    print("[FINISH]")

if __name__ == "__main__":
    main()
