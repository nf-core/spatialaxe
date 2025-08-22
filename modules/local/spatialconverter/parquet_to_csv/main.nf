process PARQUET_TO_CSV {
    tag "$meta.id"
    label 'process_low'

    container "ghcr.io/scverse/spatialdata:spatialdata0.3.0_spatialdata-io0.1.7_spatialdata-plot0.2.9"

    input:
    tuple val(meta), path(transcripts)
    val(extension)

    output:
    tuple val(meta), path("*.csv*"), emit: transcripts_csv
    path("versions.yml")           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PARQUET_TO_CSV module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    template 'parquet_to_csv.py'

    stub:
    """
    touch ${transcripts}.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatialconverter: "${task.version}"
    END_VERSIONS
    """
}
