//
// Run xeniumranger resegment
//

include { XENIUMRANGER_RESEGMENT           } from '../../modules/nf-core/xeniumranger/resegment/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow XENIUMRANGER_RESEGMENT_MORPHOLOGY {

    take:

    ch_bundle  // channel: [ val(meta), ["xenium-bundle"] ]

    main:

    ch_versions = Channel.empty()

    // run resegment with changed config values
    XENIUMRANGER_RESEGMENT ( ch_bundle )
    ch_versions = ch_versions.mix( XENIUMRANGER_RESEGMENT.out.versions )


    // run import segmentation to redine
    cells = ch_bundle.map {
        _meta, bundle -> return [ bundle + "/cells.zarr.zip" ]
    }

    XENIUMRANGER_IMPORT_SEGMENTATION (
        XENIUMRANGER_RESEGMENT.out.bundle,
        [],
        cells,
        cells,
        [],
        [],
        "pixel"
    )
    ch_versions = ch_versions.mix( XENIUMRANGER_IMPORT_SEGMENTATION.out.versions )

    emit:

    redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.bundle // channel: [ val(meta), ["redefined-xenium-bundle"] ]

    versions = ch_versions                                         // channel: [ versions.yml ]
}

