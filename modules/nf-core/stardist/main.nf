process STARDIST {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/schapirolabor/stardist:0.9.1"

    input:
    tuple val(meta), path(image)
    val model
    val maskname

    output:
    tuple val(meta), path("${prefix}/*masks.tif"), emit: mask
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def stardist_model = model ?: '2D_versatile_fluo'
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export OMP_NUM_THREADS=${task.cpus}

    # Extract DAPI channel (channel 0) from multi-channel OME-TIFF.
    # Xenium morphology_focus.ome.tif has multiple channels;
    # StarDist 2D_versatile_fluo expects single-channel input.
    python3 - ${image} <<'PYEOF'
import sys, tifffile, numpy as np
img = tifffile.imread(sys.argv[1])
orig_shape = img.shape
if img.ndim == 3:
    img = img[0]
tifffile.imwrite('dapi_input.tif', img)
print(f'Input shape: {orig_shape} -> extracted DAPI: {img.shape}')
PYEOF

    stardist-predict2d \\
        -i dapi_input.tif \\
        -o . \\
        -m ${stardist_model} \\
        ${args}

    # Convert StarDist int32 labels to uint32 for XeniumRanger compatibility
    python3 - *.stardist.tif <<'PYEOF'
import sys, tifffile, numpy as np
mask = tifffile.imread(sys.argv[1])
tifffile.imwrite(sys.argv[1], mask.astype(np.uint32))
print(f'Converted mask dtype: {mask.dtype} -> uint32, labels: {mask.max()}')
PYEOF

    mkdir -p ${prefix}
    mv *.stardist.tif ${prefix}/morphology.ome_${maskname}_masks.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stardist: \$(python -c "import stardist; print(stardist.__version__)")
        python: \$(python --version | awk '{print \$2}')
        tensorflow: \$(python -c "import tensorflow; print(tensorflow.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}
    touch ${prefix}/morphology.ome_${maskname}_masks.tif

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stardist: "0.9.1"
        python: "3.9.0"
        tensorflow: "2.10.0"
    END_VERSIONS
    """
}
