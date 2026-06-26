process SPLIT_TRANSCRIPTS {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b900c562dadb26dedce5254f88ae85440d7a08cd5e7f72cc4c3ce5aef89b5aa8/data' :
        'community.wave.seqera.io/library/pip_pandas:257725bfe0d2df83' }"

    input:
    tuple val(meta), path(transcripts)
    val(x_bins)
    val(y_bins)

    output:
    tuple val(meta), path("${meta.id}/splits.csv"), emit: splits_csv
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPLIT_TRANSCRIPTS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    utility_split_transcripts.py \\
        --transcripts ${transcripts} \\
        --x-bins ${x_bins} \\
        --y-bins ${y_bins} \\
        --prefix ${prefix}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPLIT_TRANSCRIPTS module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch "${prefix}/splits.csv"
    """
}
