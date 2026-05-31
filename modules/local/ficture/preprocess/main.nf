process FICTURE_PREPROCESS {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/08f94799a8abd47d274654c49ed5ae225811b8a64bc9788739f4c5d23fa08230/data' :
        'community.wave.seqera.io/library/pip_ficture:ad8a1265a51b53cf' }"

    input:
    tuple val(meta), path(transcripts)
    path(features)

    output:
    tuple val(meta), path("*processed_transcripts.tsv.gz"), emit: transcripts
    path("*coordinate_minmax.tsv")                        , emit: coordinate_minmax
    path("*feature.clean.tsv.gz")                         , optional:true, emit: features
    tuple val("${task.process}"), val('ficture'), eval("pip show ficture | sed -n 's/^Version: //p'"), topic: versions, emit: versions_ficture

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def features_arg = features ? "--features ${features}" : ""

    """
    ficture_preprocess.py \\
        --transcripts ${transcripts} \\
        ${features_arg} \\
        ${args}
    """

    stub:
    """
    touch processed_transcripts.tsv.gz
    touch coordinate_minmax.tsv
    """
}
