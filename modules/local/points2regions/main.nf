process POINTS2REGIONS_CLUSTER {
    tag "$meta.id"
    label 'points_cluster'
    container "community.wave.seqera.io/library/pip_points2regions:9f5bb888586554a6"

    input:
    tuple val(meta), path(transcripts)
    val(smoothing)
    val(num_clusters)

    output:
    tuple val(meta), path("clustered_s${smoothing}.csv"), emit: clustered
    tuple val(meta), path("cluster_plot_s${smoothing}.png"), emit: clustered_plot
    path "versions.yml"

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'points2regions_cluster.py'
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        points2regions_cluster: v.1.0.0
    END_VERSIONS
    """
}
