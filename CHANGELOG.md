# nf-core/spatialaxe: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.0.1dev - [date]

Initial release of nf-core/spatialaxe, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Image QC and transcript QC subworkflow (`QC`): runs `IMAGE_QC_ANALYSIS` (focus / SNR / morphology metrics) and `TRANSCRIPT_QC_PROCESSING` (per-transcript and per-cell metrics), each rendering a Quarto HTML report via the shared `QUARTO` module. New `image_qc`, `transcript_qc`, and `quarto` local modules with pinned `environment.yml` and Dockerfiles.
- GPU-optional image QC, and new QC parameters (image + transcript QC) with `nextflow_schema.json` entries.

### `Fixed`

### `Dependencies`

### `Deprecated`

## 1.0.1 - [06.08.2026]

Hotfix to tackle some bugs

### `Added`

- Template update for nf-core/tools version 4.0.3
- Adding new conf/tests folder
- Adding new test for coordinate mode to check the bugfixes
- Remove default outdir='results' for test profiles (tests)

### `Fixed`

- Only pass `--expansion-distance` to `xeniumranger import-segmentation` for nuclei-based imports. It was applied to every import, but xeniumranger rejects it for transcript-assignment imports (proseg/baysor/segger) and cells-only imports with `ERROR: --expansion-distance requires --nuclei`.
- Preserve the URI scheme (e.g. `s3://`) when building Xenium bundle child paths, so bundle validation works when the work directory is on object storage (S3/GCS/Azure). Previously `Path.toString()` dropped the scheme and the resulting path was resolved on the local filesystem, causing `Xenium bundle does not exist` / `NoSuchFileException` failures on AWS Batch.

### `Dependencies`

### `Deprecated`

## 1.0.0 - [18.06.2026]

Initial release of nf-core/spatialaxe, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
