#!/usr/bin/env python

"""Provide functions to write spatialdata object from segmentation format."""

import sys
import spatialdata
from spatialdata_io import xenium

def main():
    print("[START]")

    input_path = "${bundle}"
    output_path = "."
    outputfolder = "${outputfolder}"
    segmented_object = "${segmented_object}"

    cells_as_circles=True
    cells_boundaries=False
    nucleus_boundaries=False
    cells_labels=False
    nucleus_labels=False

    if ( segmented_object == 'cells' ):
        cells_boundaries=True
        cells_labels=True
    elif ( segmented_object == 'nuclei' ):
        nucleus_boundaries=True
        nucleus_labels=True
    elif ( segmented_object == 'cells_and_nuclei' ):
        cells_boundaries=True
        nucleus_boundaries=True
        cells_labels=True
        nucleus_labels=True
    else:
        cells_as_circles=False

    format = "${params.format}"
    if ( format == "xenium" ):
        sd_xenium_obj = xenium(
            input_path,
            cells_as_circles=cells_as_circles,
            cells_boundaries=cells_boundaries,
            nucleus_boundaries=nucleus_boundaries,
            cells_labels=cells_labels,
            nucleus_labels=nucleus_labels,
            transcripts=True,
            morphology_mip=True,
            morphology_focus=True,
            n_jobs="${task.cpus}",
        )
        print(sd_xenium_obj)
        sd_xenium_obj.write(f"{output_path}/{outputfolder}")
    else:
        sys.exit("[ERROR] Format not found")

    #Output version information
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'spatialdata: "{spatialdata.__version__}"\\n')

    print("[FINISH]")

if __name__ == "__main__":
    main()
