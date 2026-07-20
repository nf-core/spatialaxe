# nf-core/spatialaxe: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.0.0dev - [date]

Initial release of nf-core/spatialaxe, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

- Only pass `--expansion-distance` to `xeniumranger import-segmentation` for nuclei-based imports. It was applied to every import, but xeniumranger rejects it for transcript-assignment imports (proseg/baysor/segger) and cells-only imports with `ERROR: --expansion-distance requires --nuclei`.

### `Dependencies`

### `Deprecated`

## 1.0.0 - [18.06.2026]

Initial release of nf-core/spatialaxe, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
