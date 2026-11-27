process PROSEG2BAYSOR {
    tag "$meta.id"
    label 'process_high'

    container "ghcr.io/dcjones/proseg:v3.1.0"

    input:
    tuple val(meta), path(zarr_dir)

    output:
    tuple val(meta), path("${prefix}/cell-polygons.geojson")  , emit: xr_polygons
    tuple val(meta), path("${prefix}/transcript-metadata.csv"), emit: xr_metadata
    tuple val("${task.process}"), val('proseg'), eval("proseg --version | sed 's/proseg //'"), topic: versions, emit: versions_proseg

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PROSEG2BAYSOR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}

    proseg-to-baysor  \\
        ${zarr_dir} \\
        --output-transcript-metadata ${prefix}/transcript-metadata.csv \\
        --output-cell-polygons ${prefix}/cell-polygons.geojson \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PROSEG2BAYSOR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch "${prefix}/transcript-metadata.csv"
    touch "${prefix}/cell-polygons.geojson"
    """
}
