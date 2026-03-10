process CELLPOSE {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb670191b7ae1a9fd5449746453916c7014b9ea622942ca76a7cb40da7deee46/data' :
        'community.wave.seqera.io/library/python_pip_cellpose:fdf7a8c3a305a26e' }"

    input:
    tuple val(meta), path(image)
    path(model)

    output:
    tuple val(meta), path("${prefix}/*_cp_masks.tif"), emit: mask
    tuple val("${task.process}"), val('cellpose'), eval("cellpose --version | sed -n 's/cellpose version:[[:space:]]*//p' | tr -d '[:space:]'"), topic: versions, emit: versions_cellpose
    tuple val("${task.process}"), val('python'), eval("cellpose --version | sed -n 's/python version:[[:space:]]*//p' | tr -d '[:space:]'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('torch'), eval("cellpose --version | sed -n 's/torch version:[[:space:]]*//p' | tr -d '[:space:]'"), topic: versions, emit: versions_torch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def model_command = model ? "--pretrained_model ${model}" : ""
    def gpu_flag = task.accelerator ? "--use_gpu" : ""
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export OMP_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}
    # Container runs as root with HOME=/ which is not writable
    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD/.matplotlib
    export CELLPOSE_LOCAL_MODELS_PATH=\$PWD/.cellpose
    mkdir -p \$MPLCONFIGDIR \$CELLPOSE_LOCAL_MODELS_PATH

    cellpose \\
        --image_path ${image} \\
        --save_tif \\
        --verbose \\
        ${gpu_flag} \\
        ${model_command} \\
        ${args}

    # Fail fast if cellpose detected zero cells
    if grep -q "No cell pixels found" .cellpose/run.log 2>/dev/null; then
        echo "ERROR: cellpose detected 0 cells" >&2; exit 1
    fi

    mkdir -p ${prefix}
    mv *_cp_masks.tif ${prefix}/
    """

    stub:
    def name = image.name
    def base = name.lastIndexOf('.') != -1 ? name[0..name.lastIndexOf('.') - 1] : name
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}/${base}_cp_masks.tif
    """
}
