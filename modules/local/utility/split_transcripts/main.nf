process SPLIT_TRANSCRIPTS {
    tag "$meta.id"
    label 'process_low'

    container "ghcr.io/scverse/spatialdata:spatialdata0.3.0_spatialdata-io0.1.7_spatialdata-plot0.2.9"

    input:
    tuple val(meta), path(transcripts)
    val(x_bins)
    val(y_bins)

    output:
    tuple val(meta), path("splits.csv"), emit: splits_csv
    path("versions.yml")               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPLIT_TRANSCRIPTS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    template 'split_transcripts.py'

    stub:
    """
    touch ${transcripts}.parquet
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        baysor_split_parquet: "1.0.0"
    END_VERSIONS
    """
}
