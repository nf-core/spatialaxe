/*
 * DOWNSCALE_MORPHOLOGY: Pre-downscale morphology image for cellpose
 *
 * Reduces image dimensions by a scale factor so that cellpose's internal
 * rescaling (diam_mean/diameter) does not exceed GPU/CPU memory.
 * The scale factor defaults to diameter/diam_mean (e.g., 9/30 = 0.3).
 * After downscaling, cellpose should use --diameter 30 (no internal rescale).
 *
 * Input:
 *   - meta: Sample metadata map
 *   - image: Morphology OME-TIFF
 *
 * Output:
 *   - downscaled: Downscaled TIFF image
 *   - scale_info: JSON with scale factor and original dimensions
 *   - versions: Software versions
 */
process DOWNSCALE_MORPHOLOGY {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cb/cb670191b7ae1a9fd5449746453916c7014b9ea622942ca76a7cb40da7deee46/data' :
        'community.wave.seqera.io/library/python_pip_cellpose:fdf7a8c3a305a26e' }"

    input:
    tuple val(meta), path(image)

    output:
    tuple val(meta), path("${prefix}/downscaled.tif"), emit: downscaled
    tuple val(meta), path("${prefix}/scale_info.json"), emit: scale_info
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('tifffile'), eval("pip show tifffile 2>/dev/null | sed -n 's/^Version: //p'"), topic: versions, emit: versions_tifffile
    tuple val("${task.process}"), val('scikit-image'), eval("pip show scikit-image 2>/dev/null | sed -n 's/^Version: //p'"), topic: versions, emit: versions_skimage

    when:
    task.ext.when == null || task.ext.when

    script:
    def diameter = task.ext.diameter ?: 9
    def diam_mean = 30
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    utility_downscale_morphology.py \\
        --image ${image} \\
        --diameter ${diameter} \\
        --diam-mean ${diam_mean} \\
        --prefix ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/downscaled.tif
    echo '{"scale": 0.3}' > ${prefix}/scale_info.json
    """
}
