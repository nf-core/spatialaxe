process EMIT_SEGMENTATION_QC {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/94/9409ce399922a5746bea1b7df5668c3d1d79b9af49a15950d9818c4fe45ac749/data' :
        'community.wave.seqera.io/library/pandas_procs_pyarrow:d8f882b65dfea451' }"

    input:
    tuple val(meta), path(segmentation_csv)

    output:
    tuple val(meta), path("${meta.id}_segmentation_qc.json"), emit: qc_json

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    emit_segmentation_qc.py \\
        --segmentation ${segmentation_csv} \\
        --sample-id ${meta.id} \\
        --output ${meta.id}_segmentation_qc.json
    """

    stub:
    """
    echo '{"sample_id":"${meta.id}","stage":"segmentation","metrics":{"num_cells":100,"assignment_rate":0.75}}' > ${meta.id}_segmentation_qc.json
    """
}
