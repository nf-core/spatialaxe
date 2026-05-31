/*
 * EXTRACT_DAPI: Extract DAPI channel (channel 0) from multi-channel OME-TIFF.
 *
 * Xenium morphology_focus.ome.tif has multiple channels (DAPI, boundary, interior);
 * StarDist 2D_versatile_fluo expects single-channel input.
 *
 * Input:
 *   - meta: Sample metadata map
 *   - image: Multi-channel OME-TIFF morphology image
 *
 * Output:
 *   - dapi: Single-channel DAPI TIFF
 *   - versions: Software versions
 */
process EXTRACT_DAPI {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d964e0bef867bb2ff1a309c9c087d8d83ac734ce3aa315dd8311d4c1bfdafd8e/data' :
        'community.wave.seqera.io/library/python_pip_imagecodecs_nvidia-cublas-cu12_pruned:b668bcb6d531d350' }"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("${prefix}_dapi.tif"), emit: dapi
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('tifffile'), eval("python3 -c 'import tifffile; print(tifffile.__version__)'"), topic: versions, emit: versions_tifffile
    tuple val("${task.process}"), val('numpy'), eval("python3 -c 'import numpy; print(numpy.__version__)'"), topic: versions, emit: versions_numpy

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def channel_index = task.ext.channel_index ?: 0
    """
    utility_extract_dapi.py \\
        --input ${image} \\
        --output ${prefix}_dapi.tif \\
        --channel-index ${channel_index}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_dapi.tif
    """
}
