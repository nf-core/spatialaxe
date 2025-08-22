process SEGGER2XR {
    tag "$meta.id"
    label 'process_low'

    container "ghcr.io/scverse/spatialdata:spatialdata0.3.0_spatialdata-io0.1.7_spatialdata-plot0.2.9"

    input:
    tuple val(meta), path(transcripts)

    output:
    tuple val(meta), path("transcripts.parquet"), emit: transcripts_parquet
    path("versions.yml")                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SEGGER2XR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    template 'segger2xr.py'

    stub:
    """
    touch ${transcripts}.parquet
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segger2xr: "${task.version}"
    END_VERSIONS
    """
}
