# nf-core/spatialxe: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- Mode specific output:
  - [image mode](#image-mode)
  - [cooridnate mode](#coordinate-mode)
  - [segfree mode](#segfree-mode)
  - [qc mode](#qc-mode) (or using `--run_qc`)
  - [preview mode](#preview-mode)
- [Additional functionality of spatialxe](#additional-functionality):
  - [SpatialData](#spatialdata)
  - [Xenium Ranger import segmentation](#xenium-ranger-import-segmentation)
  - [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
  - [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Image mode

<details markdown="1">
<summary>Output files</summary>

- `image/`
  - `xeniumranger/`
    - `resegment/`
      - `${meta.id}/` Directory containing the output xenium bundle of Xenium
  - `baysor/`
    - `preprocess/`
      - `*.csv` filtered transcripts CSV (for Baysor 0.7.1 Parquet.jl compatibility)
    - `run/`
      - `*segmentation.csv` results of segmentation
      - `*.json` file with outlines of segmentation
      - `segmentation_params.dump.toml` file with full list of parameters used for the model
      - `segmentation_log.log` output file with metadata of running the workflow
      - `segmentation_counts.loom` loom file with metadata
      - `segmentation_cell_stats.csv` statistics of segmented cells
  - `cellpose_cells/`
    - `*masks.tif` labelled mask output from cellpose in tif format
    - `*flows.tif` cell flow output from cellpose
    - `*seg.npy` numpy array with cell segmentation data
  - `stardist_nuclei/`
    - `*.{tiff,tif}` labelled mask output from stardist in tif format
  - `resolift/`
    - `*.tiff` path to save the upscaled TIFF file

</details>

## Coordinate mode

<details markdown="1">
<summary>Output files</summary>

- `coordinate/`
  - `xenium_patch/`
    - `patches/patch_grid.json` patch_grid.json metadata file
    - `patches/patch_*/transcripts.parquet` per-patch transcripts.parquet files (one per patch)
    - `output/xr-cell-polygons.geojson` stitched cell polygons
    - `output/xr-transcript-metadata.csv` transcript metadata
  - `proseg/`
    - `preset/`
      - `cell-polygons.geojson.gz` 2D polygons for each cell in GeoJSON format. These are flattened from 3D
      - `expected-counts.csv.gz` cell-by-gene count matrix
      - `cell-metadata.csv.gz` cell centroids, volume, and other information
      - `transcript-metadata.csv.gz` transcript ids, genes, revised positions, assignment probability
      - `gene-metadata.csv.gz` per-gene summary statistics
      - `rates.csv.gz` cell-by-gene Poisson rate parameters
      - `cell-polygons-layers.geojson.gz` a separate, non-overlapping cell polygon for each z-layer, preserving 3D segmentation
      - `cell-hulls.geojson.gz` convex hulls around assigned transcripts
    - `proseg2baysor/`
      - `xr-cell-polygons.geojson` 2D polygons for each cell in GeoJSON format. These are flattened from 3D
      - `xr-transcript-metadata.csv` transcript ids, genes, revised positions, assignment probability
  - `segger/`
    - `create_dataset/`
      - `${meta.id}/` directory to save the processed Segger dataset (in PyTorch Geometric format)
    - `train/`
      - `${meta.id}/` directory to save the trained model and checkpoints
    - `predict/`
      - `${meta.id}/` directory to save the segmentation results, including cell boundaries and associations
  - `baysor/`
    - `run/`
      - `*segmentation.csv` results of segmentation
      - `*.json` file with outlines of segmentation
      - `segmentation_params.dump.toml` file with full list of parameters used for the model
      - `segmentation_log.log` output file with metadata of running the workflow
      - `segmentation_counts.loom` loom file with metadata
      - `segmentation_cell_stats.csv` statistics of segmented cells

</details>

## Segfree mode

<details markdown="1">
<summary>Output files</summary>

- `segfree/`
  - `baysor/`
    - `preprocess/`
      - `*.csv` filtered transcripts CSV (for Baysor 0.7.1 Parquet.jl compatibility)
    - `segfree/`
      - `ncvs.loom` loom file with neighborhood results
      - `ncvs_segfree_log.log` Log file with summary statistics
  - `ficture/`
    - `preprocess/`
      - `processed_transcripts.tsv.gz` transcirpt file used for FICTURE
      - `coordinate_minmax.tsv` listing the min and max of the coordinates used for FICTURE
      - `feature.clean.tsv.gz` another file contains the (unique) names of genes that should be used for FICUTRE
    - `${meta.id}/results/` files containing the results of FICTURE

</details>

## QC mode

<details markdown="1">
<summary>Output files</summary>

- `opt/`
  - `flip/`
    - `*.fa` the forward oriented fasta file
  - `track/`
    - `*.tsv` TSV file containing the gene and transcript information to which each probe aligns
  - `stat/`
    - `*.tsv` TSV file containing the summary stats
- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

## Preview mode

<details markdown="1">
<summary>Output files</summary>

- `preview/`
  - `baysor/`
    - `preview/`
      - `preview.html` segmentation preview

</details>

## Additional Functionality

### SpatialData

The pipeline create spatialdata objects (data bundles) on various stages (see metromap in the [README](../README.md))

<details markdown="1">
<summary>Output files</summary>

- `spatialdata/`
  - `write/${meta.id}/spatialdata/` spatialdata bundle of the raw data
  - `meta/${meta.id}/spatialdata_spatialxe_final/` spatialdata bundle of the final data with metadata
    - `sdata['raw_table'].uns['spatialdata_attrs']` provenance metadata
    - `sdata['raw_table'].uns['experiment_xenium']` experimental metadata
    - `sdata['raw_table'].uns['gene_panel']` gene panel metadata

</details>

### Xenium Ranger Import Segmentation)

This step is needed to import segemntations from different methods into the xenium bundle and is called at different stages of the pipeline.

<details markdown="1">
<summary>Output files</summary>

- `xeniumranger/`
  - `import_segementation/`
    - `${meta.id}/` directory containing the output xenium bundle of Xenium

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
