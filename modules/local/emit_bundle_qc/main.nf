process EMIT_BUNDLE_QC {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/94/9409ce399922a5746bea1b7df5668c3d1d79b9af49a15950d9818c4fe45ac749/data' :
        'community.wave.seqera.io/library/pandas_procs_pyarrow:d8f882b65dfea451' }"

    input:
    tuple val(meta), path(bundle)

    output:
    tuple val(meta), path("${meta.id}_${task.ext.stage ?: 'bundle_qc'}_qc.json"), emit: qc_json

    when:
    task.ext.when == null || task.ext.when

    script:
    def stage = task.ext.stage ?: 'bundle_qc'
    """
    emit_bundle_qc.py \\
        --bundle ${bundle} \\
        --sample-id ${meta.id} \\
        --stage ${stage} \\
        --output ${meta.id}_${stage}_qc.json
    """

    stub:
    def stage = task.ext.stage ?: 'bundle_qc'
    """
    echo '{"sample_id":"${meta.id}","stage":"${stage}","metrics":{"total_transcripts":50000,"xenium_num_cells":500}}' > ${meta.id}_${stage}_qc.json
    """
}
