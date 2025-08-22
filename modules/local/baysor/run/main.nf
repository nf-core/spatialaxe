process BAYSOR_RUN {
    tag "$meta.id"
    label 'process_high'

    container "khersameesh24/baysor:0.7.1"

    input:
    tuple val(meta), path(transcripts)
    path(prior_segmentation)
    path(config)
    val(scale)

    output:
    tuple val(meta), path("segmentation.csv"), emit: segmentation
    path("segmentation_polygons_2d.json")    , emit: polygons2d
    path("segmentation_polygons_3d.json")    , emit: polygons3d
    path("*.toml")                           , emit: params
    path("*.log")                            , emit: log
    path("*.loom")                           , emit: loom
    path("*.html")                           , emit: htmls
    path("segmentation_cell_stats.csv")      , emit: stats
    path("versions.yml")                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "BAYSOR_RUN module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prior_seg = "${prior_segmentation}" ? "${prior_segmentation}" : ""
    def scaling_factor = scale ? "--scale=${scale}": ""

    """
    baysor run \\
    ${transcripts} \\
    ${prior_seg} \\
    ${scaling_factor} \\
    --config=${config} \\
    --plot \\
    --polygon-format=GeometryCollectionLegacy \\
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        baysor: 0.7.1
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "BAYSOR_RUN module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    """
    touch segmentation.csv
    touch segmentation_polygons_2d.json
    touch segmentation_polygons_3d.json
    touch segmentation_log.log
    touch segmentation_counts.loom
    touch segmentation_cell_stats.csv
    touch segmentation_params.dump.toml
    touch segmentation_run.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        baysor: 0.7.1
    END_VERSIONS
    """
}
