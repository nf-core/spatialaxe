process FICTURE {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/08f94799a8abd47d274654c49ed5ae225811b8a64bc9788739f4c5d23fa08230/data' :
        'community.wave.seqera.io/library/pip_ficture:ad8a1265a51b53cf' }"

    input:
    tuple val(meta), path(transcripts)
    path(coordinate_minmax)
    path(features)

    output:
    tuple val(meta), path("results/**"), emit: results
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def features_list = features ? "--in-feature ${features}": ""

    """
    ficture run_together \\
        --in-tsv ${transcripts} \\
        --in-minmax ${coordinate_minmax} \\
        ${features_list} \\
        --out-dir results \\
        --train-width 12,18 \\
        --n-factor 6,12 \\
        --n-jobs ${task.cpus} \\
        --plot-each-factor \\
        --all

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ficture: \$(pip show ficture | grep "^Version:" | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ficture: \$(pip show ficture | grep "^Version:" | awk '{print \$2}')
    END_VERSIONS
    """
}
