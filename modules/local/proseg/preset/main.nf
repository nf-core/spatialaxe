process PROSEG {
    tag "${meta.id}"
    label 'process_high'

    container "ghcr.io/dcjones/proseg:v3.1.0"

    input:
    tuple val(meta), path(transcripts)

    output:
    tuple val(meta), path("${prefix}/cell-polygons.geojson.gz"), path("${prefix}/transcript-metadata.csv.gz"), emit: seg_outs
    tuple val(meta), path("${prefix}/proseg-output.zarr"), emit: zarr
    tuple val("${task.process}"), val('proseg'), eval("proseg --version | sed 's/proseg //'"), topic: versions, emit: versions_proseg

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("PROSEG module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // check for platform values
    if (!(params.format in ['xenium', 'cosmx', 'merscope'])) {
        error("${params.format} is an invalid platform type. Please specify xenium, cosmx, or merscope")
    }

    """
    mkdir -p ${prefix}

    proseg \\
        --${params.format} \\
        ${transcripts} \\
        --nthreads ${task.cpus} \\
        --output-expected-counts ${prefix}/expected-counts.csv.gz \\
        --output-cell-metadata ${prefix}/cell-metadata.csv.gz \\
        --output-transcript-metadata ${prefix}/transcript-metadata.csv.gz \\
        --output-gene-metadata ${prefix}/gene-metadata.csv.gz \\
        --output-rates ${prefix}/rates.csv.gz \\
        --output-cell-polygons ${prefix}/cell-polygons.geojson.gz \\
        --output-cell-polygon-layers ${prefix}/cell-polygons-layers.geojson.gz \\
        --output-union-cell-polygons ${prefix}/union-cell-polygons.geojson.gz \\
        --output-spatialdata ${prefix}/proseg-output.zarr \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("PROSEG module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/
    touch "${prefix}/expected-counts.csv.gz"
    touch "${prefix}/cell-metadata.csv.gz"
    touch "${prefix}/transcript-metadata.csv.gz"
    touch "${prefix}/gene-metadata.csv.gz"
    touch "${prefix}/rates.csv.gz"
    touch "${prefix}/cell-polygons.geojson.gz"
    touch "${prefix}/cell-polygons-layers.geojson.gz"
    touch "${prefix}/union-cell-polygons.geojson.gz"
    mkdir -p "${prefix}/proseg-output.zarr"
    """
}
