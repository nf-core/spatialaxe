//
// Run xeniumranger resegment
//

include { XENIUMRANGER_RESEGMENT           } from '../../../modules/nf-core/xeniumranger/resegment/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow XENIUMRANGER_RESEGMENT_MORPHOLOGY_OME_TIF {

    take:

    ch_bundle_path  // channel: [ val(meta), ["path-to-xenium-bundle"] ]

    main:

    ch_versions         = Channel.empty()
    ch_redefined_bundle = Channel.empty()

    // run resegment with changed config values
    XENIUMRANGER_RESEGMENT ( ch_bundle_path )
    ch_versions = ch_versions.mix( XENIUMRANGER_RESEGMENT.out.versions )


    // run import segmentation to redine xenium bundle along with nuclear segmentation
    cells = XENIUMRANGER_RESEGMENT.out.bundle.map {
        _meta, bundle -> return [ bundle + "/cells.zarr.zip" ]
    }

    // adjust the nuclear expansion distance without altering nuclei detection
    if ( params.nucleus_segmentation_only ) {

        XENIUMRANGER_IMPORT_SEGMENTATION (
            XENIUMRANGER_RESEGMENT.out.bundle,
            [],
            cells,
            [],
            [],
            [],
            "pixels"
        )
        ch_versions = ch_versions.mix( XENIUMRANGER_IMPORT_SEGMENTATION.out.versions )

        ch_redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.bundle

    } else {

        ch_redefined_bundle = XENIUMRANGER_RESEGMENT.out.bundle
    }

    emit:

    redefined_bundle = ch_redefined_bundle // channel: [ val(meta), ["redefined-xenium-bundle"] ]

    versions = ch_versions                 // channel: [ versions.yml ]
}

