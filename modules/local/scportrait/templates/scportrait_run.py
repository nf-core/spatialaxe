#!/usr/bin/env python

import os
import matplotlib.pyplot as plt
from scportrait.pipeline.project import Project
from scportrait.pipeline.extraction.hdf5 import HDF5CellExtraction

def run_scportrait(
    sdata_path,
    project_location="scportrait_output",
    cell_id_identifier=None
):
    os.makedirs(project_location, exist_ok=True)

    s = sdata.read(sdata_path)
    image_names = list(s.images.keys())
    labels = s.labels.keys()

    cyto_key = "seg_all_cytosol"
    nucleus_key = "seg_all_nucleus"

    for image_name in image_names:
        print(f"[INFO] Processing image: {image_name}")

        use_cyto = cyto_key in labels
        use_nucleus = nucleus_key in labels

        if not use_cyto and not use_nucleus:
            print(f"[WARN] Skipping {image_name} — no segmentation masks found.")
            continue

        if use_cyto and use_nucleus:
            segmentation_mask_for_extraction = cyto_key
        elif use_cyto:
            segmentation_mask_for_extraction = cyto_key
        elif use_nucleus:
            segmentation_mask_for_extraction = nucleus_key

        sub_project_location = os.path.join(project_location, image_name)
        os.makedirs(sub_project_location, exist_ok=True)

        project = Project(
            os.path.abspath(sub_project_location),
            config_path={"segmentation_mask": segmentation_mask_for_extraction},
            overwrite=True,
            debug=False,
            segmentation_f=None,
            extraction_f=HDF5CellExtraction,
            featurization_f=None,
            selection_f=None
        )

        load_kwargs = dict(
            sdata_path=sdata_path,
            input_image_name=image_name,
            cytosol_segmentation_name=cyto_key if use_cyto else None,
            nucleus_segmentation_name=nucleus_key if use_nucleus else None,
            overwrite=True,
            keep_all=True,
            remove_duplicates=True
        )

        if cell_id_identifier is not None:
            load_kwargs["cell_id_identifier"] = cell_id_identifier

        project.load_input_from_sdata(**load_kwargs)

        project.extract()

        fig1 = project.plot_segmentation_masks()
        fig1.savefig(os.path.join(sub_project_location, f"{image_name}_masks.tif"))

        fig2 = project.plot_single_cell_images()
        fig2.savefig(os.path.join(sub_project_location, f"{image_name}_cells.png"))

    # Write version info
    with open(os.path.join(project_location, "versions.yml"), "w") as f:
        f.write(f'"SCPORTRAIT":\n  scportrait: "v1.0.0"\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdata_path", required=True, help="Path to spatialdata object (Zarr or h5ad)")
    parser.add_argument("--project_location", default="scportrait_output")
    parser.add_argument("--cell_id_identifier", default=None)
    args = parser.parse_args()

    run_scportrait_all_images(
        sdata_path=args.sdata_path,
        project_location=args.project_location,
        cell_id_identifier=args.cell_id_identifier
    )
