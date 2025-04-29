//
// Run xeniumranger import-segmentation
//

include { XENIUMRANGER_IMPORT_SEGMENTATION as IMP_SEG_COUNT_MATRIX_EXP_DISTANCE } from '../../modules/nf-core/xeniumranger/import-segmentation/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION as IMP_SEG_POLYGON_GEOJSON_INPUT     } from '../../modules/nf-core/xeniumranger/import-segmentation/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION as IMP_SEG_TRANS_MATRIX_INPUT        } from '../../modules/nf-core/xeniumranger/import-segmentation/main'


workflow XENIUMRANGER_IMPORT_SEGMENTATION_REDEFINE_BUNDLE {

    take:

    ch_bundle // channel: [ val(meta), [ "xenium-bundle" ] ]

    main:

    ch_versions = Channel.empty()
    ch_redefined_bundle = Channel.empty()

    cells = ch_bundle.map {
        _meta, bundle -> return [ bundle + "/cells.zarr.zip" ]
    }

    // scenario - 1 change nuclear expansion distance / create a nucleus-only count matrix(--expansion_distance=0)
    if ( params.expansion_distance == 0 || params.expansion_distance != 5 ){

        IMP_SEG_COUNT_MATRIX_EXP_DISTANCE (
            ch_bundle,
            [],
            cells,
            [],
            [],
            [],
            "pixel"
        )
        ch_redefined_bundle = IMP_SEG_COUNT_MATRIX_EXP_DISTANCE.out.bundle
    }

    // scenario - 2 polygon input - geojson format (from QuPath)
    if ( params.qupath_polygons ) {

        IMP_SEG_POLYGON_GEOJSON_INPUT (
            ch_bundle,
            [],
            params.qupath_polygons,
            params.qupath_polygons,
            [],
            [],
            "pixel"
        )
        ch_redefined_bundle = IMP_SEG_POLYGON_GEOJSON_INPUT.out.bundle

    } else if ( params.qupath_polygons && params.nucleus_segmentation_only ) {
        IMP_SEG_POLYGON_GEOJSON_INPUT (
            ch_bundle,
            [],
            params.qupath_polygons,
            [],
            [],
            [],
            "pixel"
        )
        ch_redefined_bundle = IMP_SEG_POLYGON_GEOJSON_INPUT.out.bundle
    }

    // scenario 3 - mask input - included in the cellpose subworkflow

    // scenario 4 - transcript assignment input - included in the baysor & proseg subworkflows

    // scenario 5 - transformation matrix input
    if ( params.qupath_polygons && params.alignment_csv ) {

        IMP_SEG_TRANS_MATRIX_INPUT (
            ch_bundle,
            params.alignment_csv,
            params.qupath_polygons,
            params.qupath_polygons,
            [],
            [],
            "pixel"
        )
        ch_redefined_bundle = IMP_SEG_TRANS_MATRIX_INPUT.out.bundle
    }


    emit:

    redefined_bundle  = ch_redefined_bundle // channel: [ val(meta), ["redefined-xenium-bundle"] ]

    versions          = ch_versions         // channel: [ versions.yml ]
}

