//
// generate spatialdata object from the spatialxe layers
//

include { SPATIALDATA_WRITE as SPATIALDATA_WRITE_RAW_BUNDLE       } from '../../../modules/local/spatialdata/write/main'
include { SPATIALDATA_WRITE as SPATIALDATA_WRITE_REDEFINED_BUNDLE } from '../../../modules/local/spatialdata/write/main'
include { SPATIALDATA_MERGE as SPATIALDATA_MERGE_RAW_REDEFINED    } from '../../../modules/local/spatialdata/merge/main'
include { SPATIALDATA_META                                        } from '../../../modules/local/spatialdata/meta/main'

workflow SPATIALDATA_WRITE_META_MERGE {

    take:
    ch_bundle_path          // channel: [ val(meta), [ "path-to-xenium-bundle" ] ]
    ch_redefined_bundle     // channel: [ val(meta), [ "redefined-xenium-bundle" ] ]

    main:

    ch_versions = Channel.empty()
    ch_segmented_object = Channel.empty()

    // check segmentation - only nuclei, cells or both cells & nuclei
    if ( params.mode == 'image') {

        if ( params.nucleus_segmentation_only && params.cell_segmentation_only ) {

            ch_segmented_object = Channel.value('cells_and_nuclei')

        }

        else if ( params.nucleus_segmentation_only ) {

            ch_segmented_object = Channel.value('nuclei')

        }

        else if ( params.cell_segmentation_only ) {

            ch_segmented_object = Channel.value('cells')

        } else {

            ch_segmented_object = Channel.value([])

        }
    }

    // set all boundaries as false - default
    if ( params.mode == 'coordinate') {

        ch_segmented_object = Channel.value([])

    }


    // write spatialdata object from the raw xenium bundle
    SPATIALDATA_WRITE_RAW_BUNDLE (
        ch_bundle_path,
        'spatialdata_raw',
        ch_segmented_object
    )
    ch_versions = ch_versions.mix ( SPATIALDATA_WRITE_RAW_BUNDLE.out.versions )


    // write spatialdata object after running IMP_SEG
    SPATIALDATA_WRITE_REDEFINED_BUNDLE (
        ch_redefined_bundle,
        'spatialdata_redefined',
        ch_segmented_object
    )
    ch_versions = ch_versions.mix ( SPATIALDATA_WRITE_REDEFINED_BUNDLE.out.versions )


    // merge raw & redefined spatialdata objects
    ch_just_redefined_bundle = SPATIALDATA_WRITE_REDEFINED_BUNDLE.out.spatialdata.map {
        _meta, bundle -> return [ bundle ]
    }
    SPATIALDATA_MERGE_RAW_REDEFINED (
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        ch_just_redefined_bundle
    )
    ch_versions = ch_versions.mix ( SPATIALDATA_MERGE_RAW_REDEFINED.out.versions )


    // write metadata with spatialdata object
    ch_just_bundle_path = ch_bundle_path.map {
        _meta, bundle -> return [ bundle ]
    }
    SPATIALDATA_META (
        SPATIALDATA_MERGE_RAW_REDEFINED.out.spatialxe_bundle,
        ch_just_bundle_path
    )
    ch_versions = ch_versions.mix ( SPATIALDATA_META.out.versions )

    emit:

    ch_sd_raw       = SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata         // channel: [ val(meta), "spatialdata_raw" ]
    ch_sd_redefined = SPATIALDATA_WRITE_REDEFINED_BUNDLE.out.spatialdata   // channel: [ val(meta), "spatialdata_redefined" ]
    ch_sd_merged    = SPATIALDATA_MERGE_RAW_REDEFINED.out.spatialxe_bundle // channel: [ val(meta), "spatialdata_spatialxe" ]
    ch_sd_meta      = SPATIALDATA_META.out.spatialxe_bundle                // channel: [ val(meta), "spatialdata_spatialxe_final" ]

    versions        = ch_versions                                          // channel: [ versions.yml ]
}
