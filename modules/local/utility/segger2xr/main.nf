process SEGGER2XR {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb8fc03fa657c164c5d83f075578bbb5d9c10f1178165f94e94f33c67efca1a1/data' :
        'community.wave.seqera.io/library/spatialdata-io_spatialdata:b264928c30680e87' }"

    input:
    tuple val(meta), path(transcripts)

    output:
    tuple val(meta), path("${meta.id}/segmentation.csv")           , emit: segmentation_csv
    tuple val(meta), path("${meta.id}/transcripts.parquet")        , emit: transcripts_parquet
    tuple val(meta), path("${meta.id}/segmentation_polygons.json") , emit: viz_polygons
    path("versions.yml")                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SEGGER2XR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    template 'segger2xr.py'

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SEGGER2XR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    echo 'transcript_id,x,y,z,gene,cell,is_noise' > "${prefix}/segmentation.csv"
    touch "${prefix}/transcripts.parquet"
    echo '{"type":"FeatureCollection","features":[]}' > "${prefix}/segmentation_polygons.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        segger2xr: "${task.version}"
    END_VERSIONS
    """
}
