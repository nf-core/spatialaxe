process SEGGER_TRAIN {
    tag "$meta.id"
    label 'process_high'

    container "khersameesh24/segger:0.1.0"

    input:
    tuple val(meta), path(dataset_dir)

    output:
    tuple val(meta), path("${meta.id}_trained_models")   , emit: trained_models
    path("versions.yml")                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SEGGER_TRAIN module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def script_path = "/workspace/segger_dev/src/segger/cli/train_model.py"

    """
    python3 ${script_path} \\
        --dataset_dir ${dataset_dir} \\
        --models_dir ${prefix}_trained_models \\
        --sample_tag ${prefix} \\
        --num_workers ${task.cpus} \\
        --batch_size ${task.batch_size} \\
        --max_epochs ${task.max_epochs} \\
        --devices ${task.devices} \\
        --accelerator ${params.segger_accelerator} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segger: 0.1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_trained_models/
    touch ${prefix}_trained_models/fakefile.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segger: 0.1.0
    END_VERSIONS
    """
}
