process SEGGER_CREATE_DATASET {
    tag "${meta.id}"
    label 'process_xl'

    container "quay.io/dongzehe/segger:1.0.14"

    input:
    tuple val(meta), path(base_dir)

    output:
    tuple val(meta), path("${prefix}/"), emit: datasetdir
    tuple val("${task.process}"), val('segger'), eval("pip show segger | sed -n 's/^Version: //p'"), topic: versions, emit: versions_segger

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGGER_CREATE_DATASET module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache
    mkdir -p \$NUMBA_CACHE_DIR

    segger_create_dataset.py \\
        --bundle-dir ${base_dir} \\
        --output-dir ${prefix} \\
        --n-workers ${task.cpus} \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGGER_CREATE_DATASET module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}/
    touch "${prefix}/fake_file.txt"
    """
}
