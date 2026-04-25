process STITCH_SCS_MASKS {
    tag "${meta.id}"
    label 'process_medium'

    container "python:3.10-slim"

    input:
    tuple val(meta), path(mask_patches)

    output:
    tuple val(meta), path("cell_masks_stitched.png"), emit: stitched_mask
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pip install --no-cache-dir numpy pillow scikit-image tifffile
    
    python - <<'PYSCRIPT'
import os
import re
import numpy as np
from PIL import Image
from pathlib import Path

def parse_patch_name(filename):
    \"\"\"Parse the patch coordinates from filename like cell_masks_0:0:0:0.png\"\"\"
    match = re.search(r'cell_masks_(\\d+):(\\d+):(\\d+):(\\d+)', filename)
    if match:
        return tuple(map(int, match.groups()))
    return None

def get_image_dimensions(mask_files):
    \"\"\"Get the full image dimensions by finding max coordinates\"\"\"
    max_x, max_y = 0, 0
    patch_size = None
    
    for f in mask_files:
        coords = parse_patch_name(os.path.basename(f))
        if coords:
            patch_x, patch_y, x_idx, y_idx = coords
            max_x = max(max_x, patch_x)
            max_y = max(max_y, patch_y)
            
            if patch_size is None:
                img = Image.open(f)
                patch_size = img.size[0]
    
    return (max_x, max_y, patch_size)

def stitch_masks(mask_files, output_file):
    \"\"\"Stitch patch masks together\"\"\"
    full_width, full_height, patch_size = get_image_dimensions(mask_files)
    
    if patch_size is None:
        raise ValueError("Could not determine patch size from mask files")
    
    img_height = full_height
    img_width = full_width
    output_array = np.zeros((img_height, img_width), dtype=np.uint32)
    
    for patch_file in sorted(mask_files):
        coords = parse_patch_name(os.path.basename(patch_file))
        if not coords:
            continue
        
        patch_x, patch_y, x_idx, y_idx = coords
        patch_img = np.array(Image.open(patch_file))
        
        y_start = y_idx * patch_size
        x_start = x_idx * patch_size
        y_end = min(y_start + patch_img.shape[0], img_height)
        x_end = min(x_start + patch_img.shape[1], img_width)
        
        patch_height = y_end - y_start
        patch_width = x_end - x_start
        
        output_array[y_start:y_end, x_start:x_end] = np.maximum(
            output_array[y_start:y_end, x_start:x_end],
            patch_img[:patch_height, :patch_width]
        )
    
    output_img = Image.fromarray(output_array.astype(np.uint16))
    output_img.save(output_file)
    print(f"Stitched mask saved to {output_file}")
    print(f"Output dimensions: {output_array.shape}")

mask_files = sorted(Path('.').glob('cell_masks_*.png'))
if not mask_files:
    raise FileNotFoundError("No cell_masks_*.png files found")

mask_files = [str(f) for f in mask_files]
print(f"Found {len(mask_files)} patch files to stitch")
for f in mask_files:
    print(f"  - {f}")

stitch_masks(mask_files, "cell_masks_stitched.png")
PYSCRIPT

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pillow: \$(pip show pillow | grep Version | cut -d' ' -f2)
        scikit-image: \$(pip show scikit-image | grep Version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    touch cell_masks_stitched.png

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: 3.10
        pillow: 9.0
        scikit-image: 0.19
    END_VERSIONS
    """
}
