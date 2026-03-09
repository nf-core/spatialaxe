//
// Runs proseg with tiling: divide transcripts -> proseg per patch -> proseg2baysor -> stitch -> xeniumranger
//

include { XENIUM_PATCH_DIVIDE              } from '../../../modules/local/xenium_patch/divide/main'
include { PROSEG                           } from '../../../modules/local/proseg/preset/main'
include { PROSEG2BAYSOR                    } from '../../../modules/local/proseg/proseg2baysor/main'
include { XENIUM_PATCH_STITCH              } from '../../../modules/local/xenium_patch/stitch/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow PROSEG_PRESET_PROSEG2BAYSOR_TILED {

    take:
    ch_bundle_path         // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    ch_transcripts_parquet // channel: [ val(meta), [ "transcripts.parquet" ] ]
    ch_morphology_image    // channel: [ val(meta), ["morphology_focus.ome.tif"] ]

    main:

    ch_versions = Channel.empty()
    ch_coordinate_space = Channel.value("microns")

    // Step 1: Divide transcripts into overlapping patches
    // XENIUM_PATCH_DIVIDE expects [meta, transcripts, image] for image dimensions
    ch_divide_input = ch_transcripts_parquet
        .join(ch_morphology_image, by: 0)

    XENIUM_PATCH_DIVIDE ( ch_divide_input )
    ch_versions = ch_versions.mix( XENIUM_PATCH_DIVIDE.out.versions )

    // Step 2: Fan out patches for parallel processing
    // transpose() emits one item per patch file: [meta, parquet_path]
    ch_patches = XENIUM_PATCH_DIVIDE.out.patch_transcripts
        .transpose()
        .map { meta, parquet_file ->
            def patch_id = parquet_file.parent.name
            def patch_meta = meta.clone()
            patch_meta.sample_id = meta.id
            patch_meta.patch_id = patch_id
            patch_meta.id = "${meta.id}_${patch_id}"
            tuple(patch_meta, parquet_file)
        }

    // Step 3: Run proseg on each patch independently
    PROSEG ( ch_patches )
    ch_versions = ch_versions.mix( PROSEG.out.versions.first() )

    // Step 4: Convert proseg output to baysor format per patch
    PROSEG2BAYSOR ( PROSEG.out.zarr )
    ch_versions = ch_versions.mix( PROSEG2BAYSOR.out.versions.first() )

    // Step 5: Gather patch results per sample and reconstruct patches directory
    ch_for_stitch = PROSEG2BAYSOR.out.xr_polygons
        .join(PROSEG2BAYSOR.out.xr_metadata, by: 0)
        .map { patch_meta, geojson, csv ->
            tuple(patch_meta.sample_id, [patch_meta.patch_id, csv, geojson])
        }
        .groupTuple(by: 0)
        .map { sample_id, patch_data ->
            def sorted = patch_data.sort { it[0] }
            def patch_ids = sorted.collect { it[0] }
            def csvs = sorted.collect { it[1] }
            def geojsons = sorted.collect { it[2] }
            tuple(sample_id, patch_ids, csvs, geojsons)
        }

    ch_stitch_input = ch_for_stitch
        .join(
            XENIUM_PATCH_DIVIDE.out.grid
                .map { meta, grid -> tuple(meta.id, grid) }
        )
        .map { sample_id, patch_ids, csvs, geojsons, grid_json ->
            def meta = [id: sample_id]
            tuple(meta, grid_json, patch_ids, csvs, geojsons)
        }

    // Step 6: Reconstruct patches directory and stitch
    RECONSTRUCT_PATCHES ( ch_stitch_input )
    XENIUM_PATCH_STITCH ( RECONSTRUCT_PATCHES.out.patches_dir )
    ch_versions = ch_versions.mix( XENIUM_PATCH_STITCH.out.versions )

    // Step 7: Run xeniumranger import-segmentation
    // spatialxe signature: meta, bundle, transcript_assignment, viz_polygons, nuclei, cells, coordinate_transform, units
    ch_xr = ch_bundle_path
        .combine(XENIUM_PATCH_STITCH.out.xr_polygons_transcript, by: 0)
        .map {
            meta, bundle, xr_cell_polygons, xr_transcript_metadata -> tuple(
                meta, bundle,
                xr_transcript_metadata,
                xr_cell_polygons,
                [], [], [],
                "microns"
            )
        }

    XENIUMRANGER_IMPORT_SEGMENTATION ( ch_xr )
    ch_versions = ch_versions.mix( XENIUMRANGER_IMPORT_SEGMENTATION.out.versions_xeniumranger )

    emit:
    coordinate_space = ch_coordinate_space                          // channel: [ "microns" ]
    redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.outs    // channel: [ val(meta), ["redefined-xenium-bundle"] ]
    versions         = ch_versions                                  // channel: [ versions.yml ]
}


/*
 * RECONSTRUCT_PATCHES: Reconstruct the patches directory structure from
 * individually staged patch files for stitch_transcripts.py.
 *
 * Input:
 *   - meta: Sample metadata map
 *   - grid_json: patch_grid.json from DIVIDE
 *   - patch_ids: List of patch IDs
 *   - csv_files: Per-patch segmentation CSV files
 *   - geojson_files: Per-patch polygon GeoJSON files
 *
 * Output:
 *   - patches_dir: Reconstructed patches directory for STITCH
 *   - versions: Software versions
 */
process RECONSTRUCT_PATCHES {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f9c8f3a2de4e2aa94500011f7d7d09276e9b6f2d79ee8737c9098fe22d4649bc/data' :
        'community.wave.seqera.io/library/sopa_procps-ng_pyarrow:c9ce8cd2ede79d72' }"

    input:
    tuple val(meta), path(grid_json), val(patch_ids), path(csv_files, stageAs: 'csv_?/*'), path(geojson_files, stageAs: 'geo_?/*')

    output:
    tuple val(meta), path("patches") , emit: patches_dir
    path("versions.yml")             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ids = patch_ids instanceof List ? patch_ids : [patch_ids]
    def csvs = csv_files instanceof List ? csv_files : [csv_files]
    def geos = geojson_files instanceof List ? geojson_files : [geojson_files]
    def reconstruct_cmds = []
    for (int i = 0; i < ids.size(); i++) {
        def pid = ids[i]
        reconstruct_cmds << "mkdir -p patches/${pid}"
        reconstruct_cmds << "cp '${csvs[i]}' patches/${pid}/segmentation.csv"
        reconstruct_cmds << "cp '${geos[i]}' patches/${pid}/segmentation_polygons.json"
    }
    def reconstruct_script = reconstruct_cmds.join('\n    ')
    """
    mkdir -p patches
    cp ${grid_json} patches/patch_grid.json

    ${reconstruct_script}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | sed 's/.*version //' | sed 's/ .*//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p patches/patch_0_0
    echo '{}' > patches/patch_grid.json
    echo 'transcript_id,x,y,z,gene,cell,is_noise' > patches/patch_0_0/segmentation.csv
    echo '{"type":"FeatureCollection","features":[]}' > patches/patch_0_0/segmentation_polygons.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: "5.2.0"
    END_VERSIONS
    """
}
