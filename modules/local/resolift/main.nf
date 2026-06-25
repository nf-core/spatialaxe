process RESOLIFT {
    tag "${meta.id}"
    label 'process_low'

    container "khersameesh24/resolift:1.0.0"

    input:
    tuple val(meta), path(morphology_tiff)

    output:
    tuple val(meta), path("${prefix}/morphology.ome.enhanced.tiff"), emit: enhanced_tiff
    tuple val("${task.process}"), val('resolift'), eval("pip show resolift | sed -n 's/^Version: //p'"), topic: versions, emit: versions_resolift

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("RESOLIFT module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}

    resolift \\
        -i ${morphology_tiff} \\
        -o ${prefix}/morphology.ome.enhanced.tiff \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("RESOLIFT module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch "${prefix}/morphology.ome.enhanced.tiff"
    """
}
