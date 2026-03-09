process SEGGER_TRAIN {
    tag "${meta.id}"
    label 'process_gpu'
    maxForks params.restrict_concurrency ? 1 : 0

    container "quay.io/dongzehe/segger:1.0.14"

    input:
    tuple val(meta), path(dataset_dir)

    output:
    tuple val(meta), path("trained_models"), emit: trained_models
    path ("versions.yml"), emit: versions

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
    // Scale GPU count with retries: 4 → 8 (capped at params.devices)
    def gpu_count = Math.min((int)Math.pow(2, task.attempt + 1), params.devices as int)
    def cuda_visible = gpu_count == 1 ? "export CUDA_VISIBLE_DEVICES=0" : ""
    def accelerator = task.accelerator ? 'gpu' : 'auto'

    """
    # Set numba cache directory to avoid caching issues in container
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache
    mkdir -p \$NUMBA_CACHE_DIR

    # GPU detection logging
    echo "=== GPU Detection (SEGGER_TRAIN) ==="
    echo "Requested devices: ${gpu_count} (attempt ${task.attempt}, max ${params.devices})"
    echo "Accelerator: ${accelerator}"
    nvidia-smi 2>/dev/null && echo "GPU available: yes" || echo "GPU available: no (nvidia-smi failed)"
    python3 -c "import torch; print(f'PyTorch CUDA available: {torch.cuda.is_available()}'); print(f'CUDA device count: {torch.cuda.device_count()}')" 2>/dev/null || echo "PyTorch CUDA check failed"
    echo "===================================="

    ${cuda_visible}
    python3 ${script_path} \\
        --dataset_dir ${dataset_dir} \\
        --models_dir trained_models \\
        --sample_tag ${prefix} \\
        --batch_size ${params.batch_size_train} \\
        --max_epochs ${params.max_epochs} \\
        --devices ${gpu_count} \\
        --num_workers ${params.segger_num_workers} \\
        --accelerator ${accelerator} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segger: 0.1.0
    END_VERSIONS
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segger: 0.1.0
    END_VERSIONS
    """
}
