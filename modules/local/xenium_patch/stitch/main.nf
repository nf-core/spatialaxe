/*
 * XENIUM_PATCH_STITCH: Stitch per-patch segmentation results into unified output.
 *
 * Uses sopa's solve_conflicts() to resolve overlapping cells at patch boundaries.
 *
 * Input:
 *   - meta: Sample metadata map
 *   - patches: Directory containing patch subdirectories and patch_grid.json
 *
 * Output:
 *   - xr_polygons_transcript: Stitched cell polygons and transcript metadata
 *   - versions: Software versions
 */
process XENIUM_PATCH_STITCH {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f9c8f3a2de4e2aa94500011f7d7d09276e9b6f2d79ee8737c9098fe22d4649bc/data' :
        'community.wave.seqera.io/library/sopa_procps-ng_pyarrow:c9ce8cd2ede79d72' }"

    input:
    tuple val(meta), path(patches)

    output:
    tuple val(meta),
        path("output/xr-cell-polygons.geojson"),
        path("output/xr-transcript-metadata.csv")  , emit: xr_polygons_transcript
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('sopa'), eval('python3 -c "import sopa; print(sopa.__version__)"'), topic: versions, emit: versions_sopa

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    xenium_patch_stitch_transcripts.py \\
        --patches ${patches} \\
        --output output \\
        ${args}

    # Post-process: ensure all GeoJSON geometries are Polygon and
    # reconcile dropped cells in the transcript CSV.
    xenium_patch_stitch_postprocess.py \\
        --geojson output/xr-cell-polygons.geojson \\
        --csv output/xr-transcript-metadata.csv
    """

    stub:
    """
    mkdir -p output
    echo '{"type":"FeatureCollection","features":[]}' > output/xr-cell-polygons.geojson
    echo 'transcript_id,x,y,z,gene,cell,is_noise' > output/xr-transcript-metadata.csv
    """
}
