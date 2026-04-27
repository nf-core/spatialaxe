//
// Run SCS segmentation and import into Xenium bundle
//

include { SCS_PREPARE_MORPHOLOGY                } from '../../../subworkflows/local/scs_prepare_morphology/main'
include { SCS_SEGMENT                          } from '../../../modules/local/scs/main'
include { STITCH_SCS_MASKS                     } from '../../../modules/local/utility/stitch_scs_masks/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION     } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow SCS_IMPORT_SEGMENTATION {
    take:
    ch_morphology_image    // channel: [ val(meta), ["path-to-morphology.ome.tif"] ]
    ch_bundle_path         // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    ch_transcripts_parquet // channel: [ val(meta), ["path-to-transcripts.parquet"] ]
    ch_coordinate_space    // channel: [ val("microns") or val("pixels") ]

    main:

    ch_versions = Channel.empty()

    // Prepare SCS inputs (morphology and transcripts)
    SCS_PREPARE_MORPHOLOGY(
        ch_morphology_image,
        ch_transcripts_parquet,
    )
    ch_versions = ch_versions.mix(SCS_PREPARE_MORPHOLOGY.out.versions)

    // Combine prepared inputs for SCS_SEGMENT
    ch_scs_in = SCS_PREPARE_MORPHOLOGY.out.scs_input_bgi_tsv
        .join(SCS_PREPARE_MORPHOLOGY.out.morphology_2d, by: 0)

    // Run SCS segmentation
    SCS_SEGMENT(ch_scs_in)
    ch_versions = ch_versions.mix(SCS_SEGMENT.out.versions)

    // Stitch patch masks into single image
    STITCH_SCS_MASKS(
        SCS_SEGMENT.out.cell_masks.map { meta, masks ->
            tuple(meta, masks.flatten())
        }
    )
    ch_versions = ch_versions.mix(STITCH_SCS_MASKS.out.versions)

    // Import SCS segmentation into Xenium bundle via xeniumranger
    // Combine bundle path with stitched masks output
    ch_bundle_masks = ch_bundle_path
        .join(STITCH_SCS_MASKS.out.segmentation_csv, by: 0)
        .join(STITCH_SCS_MASKS.out.polygons, by: 0)
        .map { meta, bundle, segmentation_csv, polygons_geojson ->
            tuple(
                meta,
                bundle,
                [],
                [],
                [],
                segmentation_csv,
                polygons_geojson,
                ch_coordinate_space.val,
            )
        }

    XENIUMRANGER_IMPORT_SEGMENTATION(ch_bundle_masks)
    ch_versions = ch_versions.mix(XENIUMRANGER_IMPORT_SEGMENTATION.out.versions)

    ch_redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.bundle

    emit:
    coordinate_space = ch_coordinate_space                         // channel: [ val("microns") or val("pixels") ]
    redefined_bundle = ch_redefined_bundle                         // channel: [ val(meta), ["xenium-bundle"] ]
    versions = ch_versions                                         // channel: [ versions.yml ]
}
