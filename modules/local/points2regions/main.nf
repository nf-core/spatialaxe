process POINTS2REGIONS_CLUSTER {
    tag "$meta.id"
    label 'points_cluster'
    container "community.wave.seqera.io/library/pip_points2regions:9f5bb888586554a6"

    input:
    tuple val(meta), path(transcripts)
    val(smoothing)
    val (num_clusters)

    output:
    tuple val(meta), path("clustered_s${smoothing}.csv"), emit: clusters

    when:
    task.ext.when == null || task.ext.when

    script:

     """
    python modules/local/points2regions/templates/points2regions_cluster.py \\
        --transcripts ${transcripts} \\
        --smoothing ${smoothing} \\
        --num_clusters ${num_clusters}
    """
}

process POINTS2REGIONS_PLOT {
    tag "$meta.id"
    label 'points_visual'
    input:
    tuple val(meta), path(clusters)
    val smoothing

    output:
    path "cluster_plot.png"
    path "versions.yml"

    script:
    """
    python modules/local/points2regions/templates/points2regions_plot.py
    """
}
