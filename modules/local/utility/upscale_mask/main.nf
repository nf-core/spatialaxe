/*
 * UPSCALE_MASK: Restore cellpose masks to original image resolution
 *
 * Uses nearest-neighbor interpolation to upscale segmentation masks
 * back to original dimensions (from scale_info.json).
 *
 * Input:
 *   - meta: Sample metadata map
 *   - mask: Cellpose mask TIFF (downscaled resolution)
 *   - scale_info: JSON with original dimensions
 *
 * Output:
 *   - upscaled_mask: Mask at original resolution
 *   - versions: Software versions
 */
process UPSCALE_MASK {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb670191b7ae1a9fd5449746453916c7014b9ea622942ca76a7cb40da7deee46/data' :
        'community.wave.seqera.io/library/python_pip_cellpose:fdf7a8c3a305a26e' }"

    input:
    tuple val(meta), path(mask), path(scale_info)

    output:
    tuple val(meta), path("${prefix}/upscaled_*.tif"), emit: upscaled_mask
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('tifffile'), eval("pip show tifffile 2>/dev/null | sed -n 's/^Version: //p'"), topic: versions, emit: versions_tifffile

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    utility_upscale_mask.py \\
        --mask ${mask} \\
        --scale-info ${scale_info} \\
        --prefix ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/upscaled_mask.tif
    """
}
