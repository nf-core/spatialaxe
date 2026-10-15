process COMPUTE_PRIOR_CONFIDENCE {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::python=3.11 conda-forge::numpy conda-forge::pandas conda-forge::tifffile"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/94/9409ce399922a5746bea1b7df5668c3d1d79b9af49a15950d9818c4fe45ac749/data' :
        'community.wave.seqera.io/library/pandas_procs_pyarrow:d8f882b65dfea451' }"

    input:
    tuple val(meta), path(transcripts), path(mask)
    val(mapping)
    val(pixel_size)

    output:
    tuple val(meta), path("${meta.id}_prior_confidence_qc.json"), emit: confidence_json
    tuple val(meta), env(CONFIDENCE_VALUE),                       emit: confidence_value

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mapping_arg = mapping ? "--mapping ${mapping}" : "--mapping linear"
    """
    pip install --quiet tifffile 2>/dev/null || true

    compute_prior_confidence.py \\
        --transcripts ${transcripts} \\
        --mask ${mask} \\
        --pixel-size ${pixel_size} \\
        --sample-id ${meta.id} \\
        ${mapping_arg} \\
        --output ${meta.id}_prior_confidence_qc.json \\
        ${args}

    CONFIDENCE_VALUE=\$(python3 -c "import json; print(json.load(open('${meta.id}_prior_confidence_qc.json'))['metrics']['prior_segmentation_confidence'])")
    """

    stub:
    """
    echo '{"sample_id":"${meta.id}","stage":"prior_confidence","metrics":{"transcript_in_cell_fraction":0.65,"prior_segmentation_confidence":0.6},"mapping":"linear"}' > ${meta.id}_prior_confidence_qc.json
    CONFIDENCE_VALUE=0.6
    """
}
