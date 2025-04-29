//
// Run cellpose on the morphology tiff
//

include { CELLPOSE                         } from '../../modules/nf-core/cellpose/main'
include { RESOLIFT                         } from '../../modules/local/resolift/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow CELLPOSE_RESOLIFT_MORPHOLOGY_OME_TIF {

    take:

    ch_image  // channel: [ val(meta), ["morphology.ome.tiff"] ]
    ch_bundle // channel: [ val(meta), ["xenium-bundle"] ]

    main:

    ch_versions = Channel.empty()

    cellpose_cells = Channel.empty()
    cellpose_masks = Channel.empty()
    cellpose_flows = Channel.empty()

    // sharpen morphology tiff if param - sharpen_tiff is true
    if ( params.sharpen_tiff ) {

        RESOLIFT ( ch_image )
        ch_versions = ch_versions.mix( RESOLIFT.out.versions )

        // run cellpose on the enhanced tiff
        CELLPOSE ( RESOLIFT.out.enhanced_tiff, params.cellpose_model )
        ch_versions = ch_versions.mix( CELLPOSE.out.versions )

    } else {

        // run cellpose on the original tiff
        CELLPOSE ( ch_image, params.cellpose_model )
        ch_versions = ch_versions.mix( CELLPOSE.out.versions )
    }

    // get cellpose segmentation results
    cellpose_cells = CELLPOSE.out.cells.map {
        _meta, cells -> return [ cells ]
    }
    cellpose_masks = CELLPOSE.out.mask.map {
        _meta, mask -> return [ mask ]
    }
    cellpose_flows = CELLPOSE.out.flows.map {
        _meta, flows -> return [ flows ]
    } // TODO cellpose flows can be used as an input to the boms method

    // run import-segmentation with cellpose results
    if ( params.nucleus_segmentation_only ) {

        XENIUMRANGER_IMPORT_SEGMENTATION (
            ch_bundle,
            [],
            cellpose_masks,
            [],
            [],
            [],
            ""
        )
    } else {

        XENIUMRANGER_IMPORT_SEGMENTATION (
            ch_bundle,
            [],
            cellpose_masks,
            cellpose_cells,
            [],
            [],
            ""
        )
    }

    emit:

    mask  = CELLPOSE.out.mask                                      // channel: [ val(meta), [ "*masks.tif" ] ]
    flows = CELLPOSE.out.flows                                     // channel: [ val(meta), [ "*flows.tif" ] ]
    cells = CELLPOSE.out.cells                                     // channel: [ val(meta), [ "*seg.npy" ] ]

    redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.bundle // channel: [ val(meta), ["redefined-xenium-bundle"] ]

    versions = ch_versions                                         // channel: [ versions.yml ]
}
