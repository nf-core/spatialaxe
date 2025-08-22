process BAYSOR_PREPROCESS_TRANSCRIPTS {
    tag "$meta.id"
    label 'process_low'

    container "ghcr.io/scverse/spatialdata:spatialdata0.3.0_spatialdata-io0.1.7_spatialdata-plot0.2.9"

    input:
    tuple val(meta), path(transcripts)
    val(min_qv)
    val(max_x)
    val(min_x)
    val(max_y)
    val(min_y)

    output:
    tuple val(meta), path("*.parquet"), emit: transcripts_parquet
    path("versions.yml")              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PREPROCESS_TRANSCRIPTS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    template 'preprocess_transcripts.py'

    stub:
    """
    touch ${transcripts}.parquet
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        baysor_preprocess_transcripts: "1.0.0"
    END_VERSIONS
    """
}
