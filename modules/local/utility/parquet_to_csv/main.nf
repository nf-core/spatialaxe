process PARQUET_TO_CSV {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/94/9409ce399922a5746bea1b7df5668c3d1d79b9af49a15950d9818c4fe45ac749/data' :
        'community.wave.seqera.io/library/pandas_procs_pyarrow:d8f882b65dfea451' }"

    input:
    tuple val(meta), path(transcripts)
    val(extension)

    output:
    tuple val(meta), path("${prefix}/*.csv*"), emit: transcripts_csv
    tuple val("${task.process}"), val('pyarrow'), eval('python3 -c "import pyarrow; print(pyarrow.__version__)"'), topic: versions, emit: versions_pyarrow

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PARQUET_TO_CSV module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    utility_parquet_to_csv.py \\
        --transcripts ${transcripts} \\
        --extension ${extension} \\
        --prefix ${prefix}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PARQUET_TO_CSV module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch "${prefix}/${transcripts}.csv"
    """
}
