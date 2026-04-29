/*
 * RECONSTRUCT_PATCHES: Reconstruct the patches directory structure from
 * individually staged patch files for stitch_transcripts.py.
 */
process RECONSTRUCT_PATCHES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(grid_json), val(patch_ids), path(csv_files, stageAs: 'csv_?/*'), path(geojson_files, stageAs: 'geo_?/*')

    output:
    tuple val(meta), path("patches") , emit: patches_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    def ids = patch_ids instanceof List ? patch_ids : [patch_ids]
    def csvs = csv_files instanceof List ? csv_files : [csv_files]
    def geos = geojson_files instanceof List ? geojson_files : [geojson_files]
    """
    mkdir -p patches
    cp ${grid_json} patches/patch_grid.json

    for i in "${!ids[@]}"; do
        pid="${ids[$i]}"
        
        mkdir -p "patches/${pid}"
        cp "${csvs[$i]}" "patches/${pid}/segmentation.csv"
        cp "${geos[$i]}" "patches/${pid}/segmentation_polygons.json"
    done
    """
}
