process SEGGER_SEGMENT {
    tag "${meta.id}"
    label 'process_xl'
    label 'process_gpu'

    container "docker.io/sailmskcc/segger:1.0.0"

    input:
    tuple val(meta), path(bundle)

    output:
    tuple val(meta), path("${prefix}/segger_segmentation.parquet")                 , emit: segmentation
    tuple val(meta), path("${prefix}/export/xenium_transcript_assignment.csv")     , emit: transcript_assignment
    tuple val(meta), path("${prefix}/export/xenium_cell_boundaries.geojson")       , emit: viz_polygons
    tuple val("${task.process}"), val('segger'), eval("python -c 'from importlib.metadata import version; print(version(\"segger\"))'"), topic: versions, emit: versions_segger

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SEGGER_SEGMENT module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    if (!task.accelerator) {
        error("SEGGER_SEGMENT needs a GPU. Run with `-profile gpu` or `--use_gpu`.")
    }

    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # cupy, numba and matplotlib cache under HOME, which is read-only in the container
    export HOME=\$PWD

    segger segment \\
        -i ${bundle} \\
        -o ${prefix} \\
        --export xeniumranger \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/export
    touch ${prefix}/segger_segmentation.parquet
    echo 'transcript_id,cell,is_noise' > ${prefix}/export/xenium_transcript_assignment.csv
    echo '{"type":"FeatureCollection","features":[]}' > ${prefix}/export/xenium_cell_boundaries.geojson
    """
}
