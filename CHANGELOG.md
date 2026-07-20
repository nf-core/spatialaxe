# nf-core/spatialaxe: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.0.0dev - [date]

Initial release of nf-core/spatialaxe, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

- Preserve the URI scheme (e.g. `s3://`) when building Xenium bundle child paths, so bundle validation works when the work directory is on object storage (S3/GCS/Azure). Previously `Path.toString()` dropped the scheme and the resulting path was resolved on the local filesystem, causing `Xenium bundle does not exist` / `NoSuchFileException` failures on AWS Batch.

### `Dependencies`

### `Deprecated`

## 1.0.0 - [18.06.2026]

Initial release of nf-core/spatialaxe, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
