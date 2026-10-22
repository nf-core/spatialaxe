process TROUTPY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0e/0ebf19bee4775bb2607b945dc3583cb8c27a82a37515a08e8a34da1f19874694/data':
        'community.wave.seqera.io/library/python_pip_python-igraph_leidenalg_troutpysegmentation-freespatial-statsviz:d929c604855328d3' }"

    input:
    tuple val(meta), path(spatialdata, stageAs: "input/*")

    output:
    tuple val(meta), path("${prefix}.zarr"), emit: spatialdata
    tuple val(meta), path("report/troutpy"), emit: report
    tuple val("${task.process}"), val('troutpy'),
        eval("python -c \"import importlib.metadata as m; print(m.version('troutpy'))\""),
        topic: versions, emit: versions_troutpy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache
    export MPLCONFIGDIR=\$PWD/.mplconfig

    troutpy_urna.py \\
        --input ${spatialdata} \\
        --output-zarr ${prefix}.zarr \\
        --plots-dir report/troutpy/plots \\
        --tables-dir report/troutpy/tables \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}.zarr report/troutpy/plots report/troutpy/tables
    touch ${prefix}.zarr/.zgroup
    touch report/troutpy/plots/.stub
    touch report/troutpy/tables/.stub
    """
}
