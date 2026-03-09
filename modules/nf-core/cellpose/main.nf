process CELLPOSE {
    tag "${meta.id}"
    label 'gpu_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb670191b7ae1a9fd5449746453916c7014b9ea622942ca76a7cb40da7deee46/data' :
        'community.wave.seqera.io/library/python_pip_cellpose:fdf7a8c3a305a26e' }"

    input:
    tuple val(meta), path(image)
    val model
    val maskname

    output:
    tuple val(meta), path("${prefix}/*masks.tif"), emit: mask
    tuple val(meta), path("${prefix}/*flows.tif"), emit: flows, optional: true
    tuple val(meta), path("${prefix}/*seg.npy"), emit: cells, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("CELLPOSE module does not support conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def model_command = model ? "--pretrained_model ${model}" : ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    export NPY_PROMOTION_STATE=legacy
    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD/.matplotlib
    export CELLPOSE_LOCAL_MODELS_PATH=\$PWD/.cellpose
    mkdir -p \$MPLCONFIGDIR \$CELLPOSE_LOCAL_MODELS_PATH

    cellpose \\
        --image_path ${image} \\
        --save_tif \\
        --verbose \\
        ${model_command} \\
        ${args}

    mkdir -p ${prefix}
    mv *masks.tif ${prefix}/morphology.ome_${maskname}_masks.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(cellpose --version | awk 'NR==2 {print \$3}')
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("CELLPOSE module does not support conda. Please use Docker / Singularity / Podman instead.")
    }
    def name = image.name
    def base = name.lastIndexOf('.') != -1 ? name[0..name.lastIndexOf('.') - 1] : name
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}/morphology.ome_${maskname}_masks.tif
    touch ${prefix}/morphology.ome_${maskname}_seg.npy
    touch ${base}_cp_masks.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellpose: \$(cellpose --version | awk 'NR==2 {print \$3}')
    END_VERSIONS
    """
}
