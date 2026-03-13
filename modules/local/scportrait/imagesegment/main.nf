process SCPORTRAIT_IMAGESEGMENT {
    tag "$meta.id"
    label 'process_high'
    maxForks params.restrict_concurrency ? 1 : 0

    container "docker.io/library/python:3.11-slim"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("${prefix}/nuclei_labels.tif"), emit: nuclei, optional: true
    tuple val(meta), path("${prefix}/cells_labels.tif"), emit: cells, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def nucleusOnly = params.nucleus_segmentation_only
    template('image_segmentation.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/nuclei_labels.tif
    touch ${prefix}/cells_labels.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scportrait: stub
    END_VERSIONS
    """
}