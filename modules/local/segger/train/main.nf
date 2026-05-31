process SEGGER_TRAIN {
    tag "${meta.id}"
    label 'process_xl'
    label 'process_gpu'

    container "quay.io/dongzehe/segger:1.0.14"

    input:
    tuple val(meta), path(dataset_dir)

    output:
    tuple val(meta), path("trained_models"), emit: trained_models
    tuple val("${task.process}"), val('segger'), eval("pip show segger | sed -n 's/^Version: //p'"), topic: versions, emit: versions_segger

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGGER_TRAIN module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args = task.ext.args ?: ''
    def script_path = "/workspace/segger_dev/src/segger/cli/train_model.py"
    prefix = task.ext.prefix ?: "${meta.id}"
    def gpu_count = 2 * task.attempt
    def cuda_visible = gpu_count == 1 ? "export CUDA_VISIBLE_DEVICES=0" : ""
    def accelerator = task.accelerator ? 'gpu' : 'auto'

    """
    # Set numba cache directory to avoid caching issues in container
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache
    mkdir -p \$NUMBA_CACHE_DIR

    # GPU detection logging
    echo "=== GPU Detection (SEGGER_TRAIN) ==="
    echo "Requested devices: ${gpu_count} (attempt ${task.attempt})"
    echo "Accelerator: ${accelerator}"
    nvidia-smi 2>/dev/null && echo "GPU available: yes" || echo "GPU available: no (nvidia-smi failed)"
    python3 -c "import torch; print(f'PyTorch CUDA available: {torch.cuda.is_available()}'); print(f'CUDA device count: {torch.cuda.device_count()}')" 2>/dev/null || echo "PyTorch CUDA check failed"
    echo "===================================="

    ${cuda_visible}
    python3 ${script_path} \\
        --dataset_dir ${dataset_dir} \\
        --models_dir trained_models \\
        --sample_tag ${prefix} \\
        --devices ${gpu_count} \\
        --accelerator ${accelerator} \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGGER_TRAIN module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p trained_models/
    touch trained_models/fakefile.txt
    """
}
