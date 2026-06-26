process SPATIALDATA_MERGE {
    tag "${meta.id}"
    label 'process_high_memory'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb8fc03fa657c164c5d83f075578bbb5d9c10f1178165f94e94f33c67efca1a1/data' :
        'community.wave.seqera.io/library/spatialdata-io_spatialdata:b264928c30680e87' }"

    input:
    tuple val(meta), path(raw_bundle, stageAs: "*"), path(redefined_bundle, stageAs: "*")
    val(outputfolder)

    output:
    tuple val(meta), path("spatialdata/${prefix}/${outputfolder}"), emit: merged_bundle
    tuple val("${task.process}"), val('spatialdata'), eval("pip show spatialdata | sed -n 's/^Version: //p'"), topic: versions, emit: versions_spatialdata

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SPATIALDATA_MERGE module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    spatialdata_merge.py \\
        --raw-bundle ${raw_bundle} \\
        --redefined-bundle ${redefined_bundle} \\
        --prefix ${prefix} \\
        --output-folder ${outputfolder}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("SPATIALDATA_MERGE module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p "spatialdata/${prefix}/${outputfolder}/"
    touch "spatialdata/${prefix}/${outputfolder}/fake_file.txt"
    """
}
