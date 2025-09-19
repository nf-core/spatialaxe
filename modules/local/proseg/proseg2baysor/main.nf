process PROSEG2BAYSOR {
    tag "$meta.id"
    label 'process_high'

    container "khersameesh24/proseg:2.0.0"

    input:
    tuple val(meta), path(cell_polygons), path(transcript_metadata)

    output:
    tuple val(meta), path("${prefix}/xr-cell-polygons.geojson")     , emit: xr_polygons
    tuple val(meta), path("${prefix}/xr-transcript-metadata.csv")   , emit: xr_metadata
    tuple val(meta), path("${prefix}")                              , emit: outdir
    path("versions.yml")                                            , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PROSEG2BAYSOR (preprocess) module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    """
    mkdir -p ${prefix}
    
    proseg-to-baysor  \\
        ${transcript_metadata} \\
        ${cell_polygons} \\
        --output-transcript-metadata ${prefix}/xr-transcript-metadata.csv \\
        --output-cell-polygons ${prefix}/xr-cell-polygons.geojson \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed 's/proseg //')
    END_VERSIONS

    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PROSEG module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}

    touch ${prefix}/xr-transcript-metadata.csv
    touch ${prefix}/xr-cell-polygons.geojson

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        proseg: \$(proseg --version | sed 's/proseg //')
    END_VERSIONS
    """
}
