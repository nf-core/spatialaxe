
process RESOLVI {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/pooch_python_scviva-tools:286616b41808715c'

    input:
    tuple val(meta), path(xenium_dir)

    output:
    tuple val(meta), path("*.resolvi.h5ad"), emit: h5ad
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: '--max-epochs 100'

    """
    resolvi.py \\
        --input "${xenium_dir}" \\
        --output "${prefix}.resolvi.h5ad" \\
        ${args}

    echo '"${task.process}":' > versions.yml
    echo "  scviva-tools: \$(python -c 'import scviva; print(scviva.__version__)')" >> versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch "${prefix}.resolvi.h5ad"
    echo '"${task.process}":' > versions.yml
    echo '  scviva-tools: stub' >> versions.yml
    """
}
