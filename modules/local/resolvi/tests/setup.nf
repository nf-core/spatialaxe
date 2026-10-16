process MAKE_XENIUM {
    container 'community.wave.seqera.io/library/pooch_python_scviva-tools:286616b41808715c'
    conda "${moduleDir}/../environment.yml"

    input:
    val meta

    output:
    tuple val(meta), path('synthetic'), emit: xenium

    script:
    template 'create_xenium.py'

    stub:
    """
    mkdir synthetic
    touch synthetic/cell_feature_matrix.h5 synthetic/cells.csv.gz
    """
}
