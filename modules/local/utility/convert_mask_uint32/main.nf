/*
 * CONVERT_MASK_UINT32: Convert segmentation mask to uint32 dtype.
 *
 * XeniumRanger import-segmentation requires uint32 masks.
 * StarDist outputs int32 labels by default.
 *
 * Input:
 *   - meta: Sample metadata map
 *   - mask: Segmentation mask TIFF (any integer dtype)
 *
 * Output:
 *   - mask: uint32 segmentation mask TIFF
 *   - versions: Software versions
 */
process CONVERT_MASK_UINT32 {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d9/d964e0bef867bb2ff1a309c9c087d8d83ac734ce3aa315dd8311d4c1bfdafd8e/data' :
        'community.wave.seqera.io/library/python_pip_imagecodecs_nvidia-cublas-cu12_pruned:b668bcb6d531d350' }"

    input:
    tuple val(meta), path(mask)

    output:
    tuple val(meta), path("${prefix}_uint32_mask.tif"), emit: mask
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('tifffile'), eval("python3 -c 'import tifffile; print(tifffile.__version__)'"), topic: versions, emit: versions_tifffile
    tuple val("${task.process}"), val('numpy'), eval("python3 -c 'import numpy; print(numpy.__version__)'"), topic: versions, emit: versions_numpy

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    utility_convert_mask_uint32.py \\
        --input ${mask} \\
        --output ${prefix}_uint32_mask.tif
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_uint32_mask.tif
    """
}
