//
// Run cellpose on the morphology tiff
//

include { RESOLIFT                         } from '../../../modules/local/resolift/main'
include { DOWNSCALE_MORPHOLOGY             } from '../../../modules/local/utility/downscale_morphology/main'
include { UPSCALE_MASK as UPSCALE_CELLS    } from '../../../modules/local/utility/upscale_mask/main'
include { UPSCALE_MASK as UPSCALE_NUCLEI   } from '../../../modules/local/utility/upscale_mask/main'
include { CELLPOSE as CELLPOSE_CELLS       } from '../../../modules/nf-core/cellpose/main'
include { CELLPOSE as CELLPOSE_NUCLEI      } from '../../../modules/nf-core/cellpose/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow CELLPOSE_RESOLIFT_MORPHOLOGY_OME_TIF {
    take:
    ch_morphology_image // channel: [ val(meta), ["path-to-morphology.ome.tiff"] ]
    ch_bundle_path      // channel: [ val(meta), ["path-to-xenium-bundle"] ]

    main:

    ch_versions = Channel.empty()
    ch_imp_seg_inputs = Channel.empty()
    ch_coordinate_space = Channel.value("pixels")

    // Use empty string when no model is provided; keep as plain string for val input
    cellpose_model = params.cellpose_model ?: ''

    // sharpen morphology tiff if param - sharpen_tiff is true
    if (params.sharpen_tiff) {

        RESOLIFT(ch_morphology_image)
        ch_versions = ch_versions.mix(RESOLIFT.out.versions)

        ch_image = RESOLIFT.out.enhanced_tiff
    }
    else {

        ch_image = ch_morphology_image
    }

    // Optional pre-downscale for large images to avoid cellpose OOM
    if (params.cellpose_downscale) {

        DOWNSCALE_MORPHOLOGY(ch_image)
        ch_versions = ch_versions.mix(DOWNSCALE_MORPHOLOGY.out.versions)

        ch_cellpose_input = DOWNSCALE_MORPHOLOGY.out.downscaled
        ch_scale_info = DOWNSCALE_MORPHOLOGY.out.scale_info
    }
    else {

        ch_cellpose_input = ch_image
        ch_scale_info = Channel.empty()
    }

    // run cellpose on morphology tiff (or downscaled version)
    CELLPOSE_CELLS(ch_cellpose_input, cellpose_model, 'cells')
    ch_versions = ch_versions.mix(CELLPOSE_CELLS.out.versions)

    CELLPOSE_NUCLEI(ch_cellpose_input, 'nuclei', 'nuclei')
    ch_versions = ch_versions.mix(CELLPOSE_NUCLEI.out.versions)

    // Upscale masks back to original resolution if downscaled
    if (params.cellpose_downscale) {

        ch_cells_for_upscale = CELLPOSE_CELLS.out.mask
            .combine(ch_scale_info, by: 0)
        UPSCALE_CELLS(ch_cells_for_upscale)
        ch_versions = ch_versions.mix(UPSCALE_CELLS.out.versions)

        ch_nuclei_for_upscale = CELLPOSE_NUCLEI.out.mask
            .combine(ch_scale_info, by: 0)
        UPSCALE_NUCLEI(ch_nuclei_for_upscale)
        ch_versions = ch_versions.mix(UPSCALE_NUCLEI.out.versions)

        ch_cells_mask = UPSCALE_CELLS.out.upscaled_mask
        ch_nuclei_mask = UPSCALE_NUCLEI.out.upscaled_mask
    }
    else {

        ch_cells_mask = CELLPOSE_CELLS.out.mask
        ch_nuclei_mask = CELLPOSE_NUCLEI.out.mask
    }

    // run import-segmentation with cellpose results
    if (params.nucleus_segmentation_only) {

        ch_imp_seg_inputs = ch_bundle_path
            .combine(ch_nuclei_mask, by: 0)
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
            .combine(ch_cells_mask, by: 0)
            .combine(ch_nuclei_mask, by: 0)
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
