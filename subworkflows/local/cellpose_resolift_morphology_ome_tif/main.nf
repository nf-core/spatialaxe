//
// Run cellpose on the morphology tiff
//

include { RESOLIFT                         } from '../../../modules/local/resolift/main'
include { DOWNSCALE_MORPHOLOGY             } from '../../../modules/local/utility/downscale_morphology/main'
include { UPSCALE_MASK as UPSCALE_CELLS    } from '../../../modules/local/utility/upscale_mask/main'
include { CELLPOSE as CELLPOSE_CELLS       } from '../../../modules/nf-core/cellpose/main'
include { EXTRACT_DAPI                     } from '../../../modules/local/utility/extract_dapi/main'
include { STARDIST as STARDIST_NUCLEI      } from '../../../modules/nf-core/stardist/main'
include { CONVERT_MASK_UINT32              } from '../../../modules/local/utility/convert_mask_uint32/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow CELLPOSE_RESOLIFT_MORPHOLOGY_OME_TIF {
    take:
    ch_morphology_image // channel: [ val(meta), ["path-to-morphology.ome.tiff"] ]
    ch_bundle_path      // channel: [ val(meta), ["path-to-xenium-bundle"] ]

    main:

    ch_versions = Channel.empty()
    ch_imp_seg_inputs = Channel.empty()
    ch_coordinate_space = Channel.value("pixels")

    // Use empty list when no model is provided; path input for official cellpose module
    cellpose_model = params.cellpose_model ? file(params.cellpose_model) : []
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

    // Optional pre-downscale for large images to avoid cellpose OOM
    // Only needed when running cellpose for cells (not nucleus_segmentation_only)
    if (params.cellpose_downscale && !params.nucleus_segmentation_only) {

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
    if (!params.nucleus_segmentation_only) {
        CELLPOSE_CELLS(ch_cellpose_input, cellpose_model)
        ch_versions = ch_versions.mix(CELLPOSE_CELLS.out.versions_cellpose)
    }

    // StarDist for nuclei — extract DAPI first, then run on original resolution
    EXTRACT_DAPI(ch_image)
    ch_versions = ch_versions.mix(EXTRACT_DAPI.out.versions)

    STARDIST_NUCLEI(EXTRACT_DAPI.out.dapi, [stardist_nuclei_model, []])
    ch_versions = ch_versions.mix(STARDIST_NUCLEI.out.versions_stardist)

    // Convert StarDist mask to uint32 for XeniumRanger compatibility
    CONVERT_MASK_UINT32(STARDIST_NUCLEI.out.mask)
    ch_versions = ch_versions.mix(CONVERT_MASK_UINT32.out.versions)

    ch_nuclei_mask = CONVERT_MASK_UINT32.out.mask

    // Upscale cellpose cells mask back to original resolution if downscaled
    // StarDist nuclei mask is already at original resolution (no upscale needed)
    if (params.cellpose_downscale) {

        if (!params.nucleus_segmentation_only) {
            ch_cells_for_upscale = CELLPOSE_CELLS.out.mask
                .combine(ch_scale_info, by: 0)
            UPSCALE_CELLS(ch_cells_for_upscale)
            ch_versions = ch_versions.mix(UPSCALE_CELLS.out.versions)
            ch_cells_mask = UPSCALE_CELLS.out.upscaled_mask
        }
    }
    else {

        if (!params.nucleus_segmentation_only) {
            ch_cells_mask = CELLPOSE_CELLS.out.mask
        }
    }

    // run import-segmentation with cellpose results
    if (params.nucleus_segmentation_only) {

        ch_imp_seg_inputs = ch_bundle_path
            .combine(ch_nuclei_mask, by: 0)
            .combine(ch_coordinate_space)
            .map { meta, bundle, nuclei_seg, coord_space ->
                tuple(
                    meta,
                    bundle,
                    [],
                    [],
                    nuclei_seg,
                    [],
                    [],
                    coord_space,
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
            .combine(ch_coordinate_space)
            .map { meta, bundle, cells_seg, nuclei_seg, coord_space ->
                tuple(
                    meta,
                    bundle,
                    [],
                    [],
                    nuclei_seg,
                    cells_seg,
                    [],
                    coord_space,
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
