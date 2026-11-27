process SEGGER_PREDICT {
    tag "${meta.id}"
    label 'process_xl'
    label 'process_gpu'

    container "quay.io/dongzehe/segger:1.0.14"

    input:
    tuple val(meta), path(segger_dataset)
    path models_dir
    path transcripts

    output:
    tuple val(meta), path("benchmarks_dir"), emit: benchmarks
    tuple val(meta), path("benchmarks_dir/*/segger_transcripts.parquet"), emit: transcripts
    tuple val("${task.process}"), val('segger'), eval("pip show segger | sed -n 's/^Version: //p'"), topic: versions, emit: versions_segger

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGGER_PREDICT module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    segger_predict.py \\
        --models-dir ${models_dir} \\
        --segger-data-dir ${segger_dataset} \\
        --transcripts-file ${transcripts} \\
        --benchmarks-dir benchmarks_dir \\
        --num-workers ${task.cpus} \\
        ${args}
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
    """
}
