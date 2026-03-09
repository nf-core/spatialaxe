process BAYSOR_PREPROCESS_TRANSCRIPTS {
    tag "${meta.id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/pandas_procs_pyarrow_pip_pruned:a01d9a7721ecb2b7"

    input:
    tuple val(meta), path(transcripts)
    val min_qv
    val max_x
    val min_x
    val max_y
    val min_y

    output:
    tuple val(meta), path("${prefix}/filtered_transcripts.csv"), emit: transcripts_csv
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("BAYSOR_PREPROCESS_TRANSCRIPTS module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    template('preprocess_transcripts.py')

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("BAYSOR_PREPROCESS_TRANSCRIPTS module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}/filtered_transcripts.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        baysor_preprocess_transcripts: "1.0.0"
    END_VERSIONS
    """
}
