//
// Runs proseg for the xenium format and proseg2baysor to generate cell ploygons
//

include { PROSEG                           } from '../../../modules/local/proseg/preset/main'
include { PROSEG2BAYSOR                    } from '../../../modules/local/proseg/proseg2baysor/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow PROSEG_PRESET_PROSEG2BAYSOR {
    take:
    ch_bundle_path         // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    ch_transcripts_parquet // channel: [ val(meta), [ "transcripts.parquet" ] ]

    main:

    ch_versions = Channel.empty()
    ch_coordinate_space = Channel.value("microns")

    // run proseg with the xenium format
    PROSEG(ch_transcripts_parquet)
    ch_versions = ch_versions.mix(PROSEG.out.versions)


    // run proseg-to-baysor on the zarr output from proseg v3
    PROSEG2BAYSOR(PROSEG.out.zarr)
    ch_versions = ch_versions.mix(PROSEG2BAYSOR.out.versions)


    // run xeniumranger import-segmentation
    ch_imp_seg_inputs = ch_bundle_path
        .combine(PROSEG2BAYSOR.out.xr_metadata, by: 0)
        .combine(PROSEG2BAYSOR.out.xr_polygons, by: 0)
        .map { meta, bundle, metadata, polygons2d ->
            tuple(
                meta,
                bundle,
                metadata,
                polygons2d,
                [],
                [],
                [],
                ch_coordinate_space.val,
            )
        }

    XENIUMRANGER_IMPORT_SEGMENTATION(
        ch_imp_seg_inputs
    )

    emit:
    coordinate_space = ch_coordinate_space // channel: [ "microns" ]
    redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.outs // channel: [ val(meta), ["redefined-xenium-bundle"] ]
    versions         = ch_versions // channel: [ versions.yml ]
}
