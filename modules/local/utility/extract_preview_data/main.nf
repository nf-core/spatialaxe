process EXTRACT_PREVIEW_DATA {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c6/c6ebf365fbfd7bdde9e1453d646f45c39eddde92df5922b9881785f347bdbc2b/data' :
        'community.wave.seqera.io/library/beautifulsoup4_pandas:a3f88f59088edad5' }"

    input:
    tuple val(meta), path(preview_html)

    output:
    tuple val(meta), path("${prefix}/*_mqc.tsv"), emit: mqc_data
    tuple val(meta), path("${prefix}/*_mqc.png"), emit: mqc_img
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("EXTRACT_PREVIEW_DATA module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    utility_extract_preview_data.py \\
        --preview-html ${preview_html} \\
        --prefix ${prefix}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("EXTRACT_PREVIEW_DATA module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}/noise_distribution_mqc.tsv
    touch ${prefix}/gene_structure_mqc.tsv
    touch ${prefix}/umap_mqc.tsv
    touch ${prefix}/transcript_plots_mqc.png
    touch ${prefix}/noise_level_mqc.png
    """
}
