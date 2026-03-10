process SEGGER_PREDICT {
    tag "${meta.id}"
    label 'process_gpu'

    container "quay.io/dongzehe/segger:1.0.14"

    input:
    tuple val(meta), path(segger_dataset)
    path models_dir
    path transcripts

    output:
    tuple val(meta), path("benchmarks_dir"), emit: benchmarks
    tuple val(meta), path("benchmarks_dir/*/segger_transcripts.parquet"), emit: transcripts
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGGER_PREDICT module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args = task.ext.args ?: ''
    def script_path = "/workspace/segger_dev/src/segger/cli/predict_fast.py"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Limit cupy GPU memory to 80% so PyTorch has headroom for graph attention ops
    export CUPY_GPU_MEMORY_LIMIT="80%"
    # Belt-and-suspenders: ensure PyTorch uses expandable segments (also set in env {} block)
    export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True,max_split_size_mb:512"

    # Set numba cache directory to avoid caching issues in container
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache
    mkdir -p \$NUMBA_CACHE_DIR

    # GPU detection logging
    echo "=== GPU Detection (SEGGER_PREDICT) ==="
    nvidia-smi 2>/dev/null && echo "GPU available: yes" || echo "GPU available: no (nvidia-smi failed)"
    python3 -c "import torch; print(f'PyTorch CUDA available: {torch.cuda.is_available()}'); print(f'CUDA device count: {torch.cuda.device_count()}')" 2>/dev/null || echo "PyTorch CUDA check failed"
    echo "======================================"

    # Use all available GPUs (autocast reduces VRAM ~50%, so multi-GPU is safe)
    GPU_IDS=\$(python3 -c "
import torch
n = torch.cuda.device_count()
print(','.join(str(i) for i in range(n)) if n > 0 else '0')
" 2>/dev/null || echo "0")
    echo "Using GPUs: \$GPU_IDS"

    # Patch predict_parquet.py at runtime (avoids Docker rebuild)
    PRED_PY=\$(python3 -c "import segger.prediction.predict_parquet as m; print(m.__file__)")

    # 1. Add torch.no_grad() to disable gradient graphs during inference (~30-50% VRAM savings)
    sed -i 's/with cp.cuda.Device(gpu_id):/with cp.cuda.Device(gpu_id), torch.no_grad():/' "\$PRED_PY"

    # 2. Seed random for deterministic GPU assignment (avoids stochastic OOM)
    sed -i 's/gpu_id = random.choice(gpu_ids)/random.seed(0); gpu_id = random.choice(gpu_ids)/' "\$PRED_PY"
    echo "Patched \$PRED_PY: torch.no_grad() + round-robin GPU assignment"

    python3 ${script_path} \\
        --models_dir ${models_dir} \\
        --segger_data_dir ${segger_dataset} \\
        --transcripts_file ${transcripts} \\
        --benchmarks_dir benchmarks_dir \\
        --batch_size ${params.batch_size_predict} \\
        --use_cc ${params.cc_analysis} \\
        --knn_method ${params.segger_knn_method} \\
        --num_workers ${task.cpus} \\
        --gpu_ids \$GPU_IDS \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segger: 0.1.0
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGGER_PREDICT module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p "benchmarks_dir"
    touch "benchmarks_dir/fake_file.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segger: 0.1.0
    END_VERSIONS
    """
}
