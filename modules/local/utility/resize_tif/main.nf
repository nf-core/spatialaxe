process RESIZE_TIF {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6d/6d5aedb8fcf066eecd9f0dfac93bfffc8161bdae65b4502509d9953db2036a7e/data' :
        'community.wave.seqera.io/library/numpy_pandas_pyarrow_scikit-image_tifffile:131397039376b375' }"

    input:
    tuple val(meta), path(transcripts), path(mask), path(metadata)

    output:
    tuple val(meta), path("${meta.id}/resized_*.tif"), emit: resized_mask
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('tifffile'), eval('python3 -c "import tifffile; print(tifffile.__version__)"'), topic: versions, emit: versions_tifffile

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("RESIZE_TIF module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    prefix = task.ext.prefix ?: "${meta.id}"

    """
    utility_resize_tif.py \\
        --mask ${mask} \\
        --transcripts ${transcripts} \\
        --metadata ${metadata} \\
        --prefix ${prefix} \\
        --mask-filename ${mask}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("RESIZE_TIF module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch "${prefix}/resized_${mask}.tif"
    """
}
