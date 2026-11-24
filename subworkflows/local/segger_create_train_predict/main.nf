//
// Run segger segment and import its segmentation with xeniumranger
//

include { SEGGER_SEGMENT                   } from '../../../modules/local/segger/segment/main'
include { XENIUMRANGER_IMPORTSEGMENTATION  } from '../../../modules/nf-core/xeniumranger/importsegmentation/main'

workflow SEGGER_CREATE_TRAIN_PREDICT {
    take:
    ch_bundle              // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    expansion_distance     // value: nuclear expansion distance

    main:

    // Xenium transcript coordinates are in microns
    ch_coordinate_space = channel.value("microns")

    SEGGER_SEGMENT(ch_bundle)

    ch_imp_seg_inputs = ch_bundle
        .join(SEGGER_SEGMENT.out.transcript_assignment)
        .join(SEGGER_SEGMENT.out.viz_polygons)
        .map { meta, bundle, transcript_assignment, polygons ->
            tuple(
                meta,
                bundle,
                transcript_assignment,
                polygons,
                [],  // nuclei
                [],  // cells
                [],  // coordinate_transform
                ch_coordinate_space.val,
                expansion_distance,
            )
        }

    XENIUMRANGER_IMPORTSEGMENTATION(
        ch_imp_seg_inputs
    )

    emit:
    coordinate_space = ch_coordinate_space // channel: [ "microns" ]
    redefined_bundle = XENIUMRANGER_IMPORTSEGMENTATION.out.outs // channel: [ val(meta), ["redefined-xenium-bundle"] ]
}
