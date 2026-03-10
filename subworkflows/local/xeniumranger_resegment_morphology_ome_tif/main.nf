//
// Run xeniumranger resegment
//

include { XENIUMRANGER_RESEGMENT           } from '../../../modules/nf-core/xeniumranger/resegment/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow XENIUMRANGER_RESEGMENT_MORPHOLOGY_OME_TIF {
    take:
    ch_bundle_path // channel: [ val(meta), ["path-to-xenium-bundle"] ]

    main:

    ch_versions = Channel.empty()
    ch_redefined_bundle = Channel.empty()
    ch_coordinate_space = Channel.value("pixels")

    // run resegment with changed config values
    XENIUMRANGER_RESEGMENT(ch_bundle_path)


    // run import segmentation to redine xenium bundle along with nuclear segmentation
    // Keep meta in the cells channel for proper per-sample joining
    def cells = XENIUMRANGER_RESEGMENT.out.outs.map { meta, bundle ->
        return [meta, bundle + "/cells.zarr.zip"]
    }

    // adjust the nuclear expansion distance without altering nuclei detection
    if (params.nucleus_segmentation_only) {

        def ch_imp_seg_inputs = ch_bundle_path
            .join(XENIUMRANGER_RESEGMENT.out.outs, by: 0)
            .join(cells, by: 0)
            .map { meta, bundle, reseg_bundle, cells_zarr ->
                tuple(
                    meta,
                    bundle,
                    [],
                    [],
                    [],
                    cells_zarr,
                    [],
                    "pixels",
                )
            }

        XENIUMRANGER_IMPORT_SEGMENTATION(
            ch_imp_seg_inputs
        )

        ch_redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.outs
    }
    else {

        ch_redefined_bundle = XENIUMRANGER_RESEGMENT.out.outs
    }

    emit:
    redefined_bundle = ch_redefined_bundle // channel: [ val(meta), ["redefined-xenium-bundle"] ]
    coordinate_space = ch_coordinate_space // channel: [ ["pixels"] ]
    versions         = ch_versions // channel: [ versions.yml ]
}
