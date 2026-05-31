/*
 * XENIUM_PATCH_DIVIDE: Split transcripts.parquet into overlapping patches.
 *
 * Input:
 *   - meta: Sample metadata map
 *   - transcripts: transcripts.parquet file
 *   - image: morphology image (for getting dimensions)
 *
 * Output:
 *   - grid: patch_grid.json metadata file
 *   - patch_transcripts: per-patch transcripts.parquet files (one per patch)
 *   - versions: Software versions
 */
process XENIUM_PATCH_DIVIDE {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f9c8f3a2de4e2aa94500011f7d7d09276e9b6f2d79ee8737c9098fe22d4649bc/data' :
        'community.wave.seqera.io/library/sopa_procps-ng_pyarrow:c9ce8cd2ede79d72' }"

    input:
    tuple val(meta), path(transcripts), path(image)

    output:
    tuple val(meta), path("patches/patch_grid.json")              , emit: grid
    tuple val(meta), path("patches/patch_*/transcripts.parquet")  , emit: patch_transcripts
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pyarrow'), eval('python3 -c "import pyarrow; print(pyarrow.__version__)"'), topic: versions, emit: versions_pyarrow

    when:
    task.ext.when == null || task.ext.when

    script:
    def tile_width = task.ext.tile_width ?: 2000
    def overlap = task.ext.overlap ?: 50
    def balanced = task.ext.balanced
    def balanced_flag = balanced == true || balanced == 'true' ? '--balanced' : ''
    """
    divide_transcripts.py \\
        --transcripts ${transcripts} \\
        --output patches \\
        --tile-width ${tile_width} \\
        --overlap ${overlap} \\
        ${balanced_flag} \\
        --image-width \$(python3 -c "import tifffile; print(tifffile.imread('${image}').shape[-1])") \\
        --image-height \$(python3 -c "import tifffile; print(tifffile.imread('${image}').shape[-2])")

    """

    stub:
    """
    mkdir -p patches/patch_0_0
    touch patches/patch_0_0/transcripts.parquet
    echo '{}' > patches/patch_grid.json
    """
}
