process SPLIT_PURIFY {
    tag "$meta.id"
    label 'process_medium'

    // Build locally from modules/local/split/Dockerfile until an image is published.
    container 'docker.io/library/split:0.3.0'

    input:
    tuple val(meta), path(counts), path(reference), path(weights)

    output:
    tuple val(meta), path("${prefix}/matrix.mtx.gz"), path("${prefix}/features.tsv.gz"), path("${prefix}/barcodes.tsv.gz"), emit: purified_counts
    tuple val(meta), path("${prefix}/cell_metadata.tsv"), emit: metadata
    tuple val("${task.process}"), val('split'), eval("Rscript -e 'cat(as.character(packageVersion(\"SPLIT\")))'"), topic: versions, emit: versions_split

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}

    split-purify \\
        --counts "${counts}" \\
        --reference "${reference}" \\
        --weights "${weights}" \\
        --outdir "${prefix}" \\
        ${args}

    gzip -n "${prefix}/matrix.mtx" "${prefix}/features.tsv" "${prefix}/barcodes.tsv"
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}"
    touch "${prefix}/cell_metadata.tsv"
    printf '' | gzip -n > "${prefix}/matrix.mtx.gz"
    printf '' | gzip -n > "${prefix}/features.tsv.gz"
    printf '' | gzip -n > "${prefix}/barcodes.tsv.gz"
    """
}
