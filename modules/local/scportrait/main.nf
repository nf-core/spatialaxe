process SCPORTRAIT {
    tag "$meta.id"
    label 'sdata'
    container "community.wave.seqera.io/library/pip_scportrait:e0651f1fbb601e73"

    input:
    tuple val(meta), path(sdata)
    val(cell_id_identifier) optional true

    output:
    tuple val(meta), path("**/*_masks.tif"), emit: sc_mask
    tuple val(meta), path("**/*_images.png"), emit: sc_images
    path "versions.yml"                 ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'scportrait_run.py'

    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scportrait_run: v.1.0.0
    END_VERSIONS
    """
}
