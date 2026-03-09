//
// Run stardist on the morphology tiff
//

include { RESOLIFT                         } from '../../../modules/local/resolift/main'
include { STARDIST as STARDIST_CELLS       } from '../../../modules/nf-core/stardist/main'
include { STARDIST as STARDIST_NUCLEI      } from '../../../modules/nf-core/stardist/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow STARDIST_RESOLIFT_MORPHOLOGY_OME_TIF {
    take:
    ch_morphology_image // channel: [ val(meta), ["path-to-morphology.ome.tiff"] ]
    ch_bundle_path      // channel: [ val(meta), ["path-to-xenium-bundle"] ]

    main:

    ch_versions = Channel.empty()
    ch_imp_seg_inputs = Channel.empty()
    ch_coordinate_space = Channel.value("pixels")

    // Use default model when no model is provided
    stardist_model = params.stardist_model ?: '2D_versatile_fluo'
    stardist_nuclei_model = params.stardist_nuclei_model ?: '2D_versatile_fluo'

    // sharpen morphology tiff if param - sharpen_tiff is true
    if (params.sharpen_tiff) {

        RESOLIFT(ch_morphology_image)
        ch_versions = ch_versions.mix(RESOLIFT.out.versions)

        ch_image = RESOLIFT.out.enhanced_tiff
    }
    else {

        ch_image = ch_morphology_image
    }

    // run stardist on morphology tiff
    STARDIST_CELLS(ch_image, stardist_model, 'cells')
    ch_versions = ch_versions.mix(STARDIST_CELLS.out.versions)

    STARDIST_NUCLEI(ch_image, stardist_nuclei_model, 'nuclei')
    ch_versions = ch_versions.mix(STARDIST_NUCLEI.out.versions)


    // run import-segmentation with stardist results
    if (params.nucleus_segmentation_only) {

        ch_imp_seg_inputs = ch_bundle_path
            .combine(STARDIST_NUCLEI.out.mask, by: 0)
            .map { meta, bundle, nuclei_seg ->
                tuple(
                    meta,
                    bundle,
                    [],
                    [],
                    nuclei_seg,
                    [],
                    [],
                    ch_coordinate_space.val,
                )
            }
        XENIUMRANGER_IMPORT_SEGMENTATION(
            ch_imp_seg_inputs
        )
    }
    else {

        ch_imp_seg_inputs = ch_bundle_path
            .combine(STARDIST_CELLS.out.mask, by: 0)
            .combine(STARDIST_NUCLEI.out.mask, by: 0)
            .map { meta, bundle, cells_seg, nuclei_seg ->
                tuple(
                    meta,
                    bundle,
                    [],
                    [],
                    nuclei_seg,
                    cells_seg,
                    [],
                    ch_coordinate_space.val,
                )
            }
        XENIUMRANGER_IMPORT_SEGMENTATION(
            ch_imp_seg_inputs
        )
    }

    emit:
    coordinate_space = ch_coordinate_space // channel: [ ["pixels"] ]
    redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.outs // channel: [ val(meta), ["redefined-xenium-bundle"] ]
    versions         = ch_versions // channel: [ versions.yml ]
}
