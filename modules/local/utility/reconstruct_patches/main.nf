/*
 * RECONSTRUCT_PATCHES: Reconstruct the patches directory structure from
 * individually staged patch files for stitch_transcripts.py.
 *
 * Inputs:
 *   meta           - sample metadata map
 *   grid_json      - patch_grid.json from XENIUM_PATCH_DIVIDE
 *   patch_ids      - list of patch identifiers (e.g. patch_0000, patch_0001, ...)
 *   csv_files      - per-patch Baysor segmentation.csv files (staged into csv_?/ dirs)
 *   geojson_files  - per-patch Baysor segmentation_polygons.json files (staged into geo_?/ dirs)
 *
 * Outputs:
 *   patches_dir    - reconstructed patches/ directory containing patch_grid.json plus
 *                    one subdirectory per patch with segmentation.csv and segmentation_polygons.json
 *   versions       - topic-channel version emission for coreutils (cp)
 */
process RECONSTRUCT_PATCHES {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b9/b900c562dadb26dedce5254f88ae85440d7a08cd5e7f72cc4c3ce5aef89b5aa8/data' :
        'community.wave.seqera.io/library/pip_pandas:257725bfe0d2df83' }"

    input:
    tuple val(meta), path(grid_json), val(patch_ids), path(csv_files, stageAs: 'csv_?/*'), path(geojson_files, stageAs: 'geo_?/*')

    output:
    tuple val(meta), path("patches"), emit: patches_dir
    tuple val("${task.process}"), val('coreutils'), eval("cp --version | head -n1 | awk '{print \$NF}'"), topic: versions, emit: versions_coreutils

    when:
    task.ext.when == null || task.ext.when

    script:
    def ids = patch_ids instanceof List ? patch_ids : [patch_ids]
    def csvs = csv_files instanceof List ? csv_files : [csv_files]
    def geos = geojson_files instanceof List ? geojson_files : [geojson_files]

    def reconstruct_script = ids.withIndex().collect { pid, idx ->
        [
            "mkdir -p patches/${pid}",
            "cp '${csvs[idx]}' patches/${pid}/segmentation.csv",
            "cp '${geos[idx]}' patches/${pid}/segmentation_polygons.json",
        ].join('\n    ')
    }.join('\n    ')
    """
    mkdir -p patches
    cp '${grid_json}' patches/patch_grid.json

    ${reconstruct_script}
    """

    stub:
    def ids = patch_ids instanceof List ? patch_ids : [patch_ids]
    def stub_files = ids.collect { pid ->
        [
            "mkdir -p patches/${pid}",
            "touch patches/${pid}/segmentation.csv",
            "touch patches/${pid}/segmentation_polygons.json",
        ].join('\n    ')
    }.join('\n    ')
    """
    mkdir -p patches
    touch patches/patch_grid.json
    ${stub_files}
    """
}
