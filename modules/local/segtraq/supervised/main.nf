process SEGTRAQ_SUPERVISED {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/priyal_tripathi/segtraq:0.0.3"

    input:
    tuple val(meta), path(spatialdata_zarr)
    path  markers
    val   cell_type_key

    output:
    tuple val(meta), path("segtraq_qc/${prefix}/"), emit: qc_results
    path("versions.yml")                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGTRAQ_SUPERVISED module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    template 'supervised.py'

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGTRAQ_SUPERVISED module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p "segtraq_qc/${prefix}"
    touch "segtraq_qc/${prefix}/supervised_summary.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segtraq: \$(pip show segtraq | grep Version | cut -d' ' -f2)
    END_VERSIONS
    """
}
