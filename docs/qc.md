# nf-core/spatialaxe: QC

## Introduction

This document describes the QC reports produced by the pipeline.

> [!NOTE]
> This is a work in progress.

## QC methods overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and inspects the quality of the data using the following methods:

- [MultiQC Xenium Extra Plugin](https://github.com/MultiQC/xenium-extra)
- [OPT](https://github.com/JEFworks-Lab/off-target-probe-tracker)
- [spoQC](https://spoqc.readthedocs.io/en/latest/index.html)
- GrandQC, when enabled with `--run_grandqc`

## MultiQC Xenium Extra Plugin

## OPT

## spoQC

A full description of the report can be found here:

- [Base report](https://spoqc.readthedocs.io/en/latest/output_report.html)
- [Additional report](https://spoqc.readthedocs.io/en/latest/output_report_additional.html)

## GrandQC

GrandQC can be enabled independently of the existing QC layer with
`--run_grandqc`. It processes each sample's morphology image after running
tissue detection, then artifact detection. The image must be a supported RGB
TIFF or OME-TIFF. OME physical pixel-size metadata is used when available;
otherwise, provide `--grandqc_mpp` in microns per pixel.

Provide the tissue-detection and artifact-detection model checkpoints with
`--grandqc_tissue_model` and `--grandqc_qc_model`. The checkpoints are not
distributed with SpatialAxe. GrandQC is released under the
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
License](../data/grandqc/LICENSE); use is subject to its terms, including its
non-commercial restriction. Verify that your use is permitted.

If you use GrandQC in research, cite: Weng Z. et al. "GrandQC: A comprehensive
solution to quality control problem in digital pathology." *Nature
Communications* (2024). https://doi.org/10.1038/s41467-024-54769-y

The module currently runs on CPU and uses a pinned Python environment to build
its runtime container through Nextflow Wave. For example:

```bash
nextflow run nf-core/spatialaxe \
    -profile docker,wave \
    --input samplesheet.csv \
    --outdir results \
    --mode qc \
    --run_grandqc \
    --grandqc_tissue_model /path/to/Tissue_Detection_MPP10.pth \
    --grandqc_qc_model /path/to/GrandQC_MPP15.pth \
    --grandqc_mpp_model 1.5
```

Add `--grandqc_mpp <microns-per-pixel>` only when the image does not contain
OME physical pixel-size metadata. The image path in the samplesheet must point
to an RGB TIFF/OME-TIFF supported by the adapted GrandQC reader.
