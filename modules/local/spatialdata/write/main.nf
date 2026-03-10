process SPATIALDATA_WRITE {
    tag "${meta.id}"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb8fc03fa657c164c5d83f075578bbb5d9c10f1178165f94e94f33c67efca1a1/data' :
        'community.wave.seqera.io/library/spatialdata-io_spatialdata:b264928c30680e87' }"

    input:
    tuple val(meta), path(bundle, stageAs: "*")
    val(outputfolder)
    val(segmented_object)
    val(coordinate_space)

    output:
    tuple val(meta), path("spatialdata/${prefix}/${outputfolder}"), emit: spatialdata
    path ("versions.yml")                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit(1, "SPATIALDATA_WRITE module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    template('write.py')

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit(1, "SPATIALDATA_WRITE module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p "spatialdata/${prefix}/${outputfolder}"
    touch "spatialdata/${prefix}/${outputfolder}/fake_file.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spatialdata: \$(echo \$( python -c "import spatialdata; print(spatialdata.__version__)" 2>&1) )
    END_VERSIONS
    """
}
