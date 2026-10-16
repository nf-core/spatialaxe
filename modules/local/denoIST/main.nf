process DENOIST {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-denoist=1.0.0"

    input:
    tuple val(meta), path(mat), path(tx), path(coords, optional: true)

    output:
    tuple val(meta), path("*_denoised.rds"), emit: denoised_data
    path "versions.yml", emit: versions

    script:
    def coords_arg = coords ? "--coords ${coords}" : ""
    """
    denoist.R \\
        --mat ${mat} \\
        --tx ${tx} \\
        ${coords_arg} \\
        --output ${meta.id}_denoised.rds \\
        --tx-x ${params.denoist_tx_x} \\
        --tx-y ${params.denoist_tx_y} \\

    printf '"${task.process}":\\n    denoist: 1.0.0\\n' > versions.yml
    """
}
