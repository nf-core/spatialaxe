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
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def diameter = task.ext.diameter ?: 9
    def diam_mean = 30
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    python3 -c "
import tifffile, numpy as np, json
from PIL import Image

diameter = ${diameter}
diam_mean = ${diam_mean}
scale = diameter / diam_mean  # e.g., 9/30 = 0.3

with tifffile.TiffFile('${image}') as tif:
    img = tif.pages[0].asarray()

orig_h, orig_w = img.shape[:2]
new_h = max(int(orig_h * scale), 256)
new_w = max(int(orig_w * scale), 256)

print(f'Original: {img.shape}, dtype={img.dtype}')
print(f'Downscaling by {scale:.3f}: ({orig_h}, {orig_w}) -> ({new_h}, {new_w})')

pil_img = Image.fromarray(img)
pil_img = pil_img.resize((new_w, new_h), Image.LANCZOS)
img_ds = np.array(pil_img, dtype=img.dtype)

tifffile.imwrite('${prefix}/downscaled.tif', img_ds, compression='zlib')
json.dump({
    'scale': scale,
    'orig_h': orig_h,
    'orig_w': orig_w,
    'new_h': new_h,
    'new_w': new_w,
    'diameter': diameter,
    'diam_mean': diam_mean
}, open('${prefix}/scale_info.json', 'w'))
print(f'Done: downscaled.tif written')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
        tifffile: \$(python3 -c 'import tifffile; print(tifffile.__version__)')
        pillow: \$(python3 -c 'import PIL; print(PIL.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/downscaled.tif
    echo '{"scale": 0.3}' > ${prefix}/scale_info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
