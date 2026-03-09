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
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}

    python3 -c "
import tifffile, numpy as np, json
from PIL import Image

info = json.load(open('${scale_info}'))
orig_h, orig_w = info['orig_h'], info['orig_w']

mask = tifffile.imread('${mask}')
print(f'Mask: {mask.shape}, dtype={mask.dtype}, unique cells: {len(np.unique(mask)) - 1}')
print(f'Upscaling to ({orig_h}, {orig_w})')

pil_mask = Image.fromarray(mask)
pil_mask = pil_mask.resize((orig_w, orig_h), Image.NEAREST)
mask_up = np.array(pil_mask, dtype=mask.dtype)

out_name = '${prefix}/upscaled_${mask.baseName}.tif'
tifffile.imwrite(out_name, mask_up, compression='zlib')
print(f'Done: {out_name}, unique cells: {len(np.unique(mask_up)) - 1}')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
        tifffile: \$(python3 -c 'import tifffile; print(tifffile.__version__)')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/upscaled_mask.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
    END_VERSIONS
    """
}
