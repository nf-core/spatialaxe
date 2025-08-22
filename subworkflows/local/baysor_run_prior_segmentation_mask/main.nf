//
// Run baysor run & import-segmentation
//

include { BAYSOR_PREPROCESS_TRANSCRIPTS    } from '../../../modules/local/baysor/preprocess/main'
include { BAYSOR_RUN                       } from '../../../modules/local/baysor/run/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'


workflow BAYSOR_RUN_PRIOR_SEGMENTATION_MASK {

    take:

    ch_bundle_path         // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    ch_transcripts_parquet // channel: [ val(meta), ["path-to-transcripts.parquet"] ]
    ch_segmentation_mask   // channel: [ ["path-to-prior-segmentation-mask"] ]
    ch_config              // channel: [ "path-to-xenium.toml" ]

    main:

    ch_versions             = Channel.empty()

    ch_transcripts          = Channel.empty()

    ch_segmentation         = Channel.empty()
    ch_polygons2d           = Channel.empty()
    ch_htmls                = Channel.empty()

    ch_redefined_bundle     = Channel.empty()

    // filter transcripts.parquet based on thresholds
    if ( params.filter_transcripts ) {

        BAYSOR_PREPROCESS_TRANSCRIPTS (
            ch_transcripts_parquet,
            params.min_qv,
            params.max_x,
            params.min_x,
            params.max_y,
            params.min_y
        )
        ch_versions = ch_versions.mix ( BAYSOR_PREPROCESS_TRANSCRIPTS.out.versions )

        ch_transcripts = BAYSOR_PREPROCESS_TRANSCRIPTS.out.transcripts_parquet

    } else {

        ch_transcripts = ch_transcripts_parquet
    }

    // run baysor with morphology.tiff
    BAYSOR_RUN (
        ch_transcripts,
        ch_segmentation_mask,
        ch_config,
        30
    )
    ch_versions = ch_versions.mix( BAYSOR_RUN.out.versions )

    ch_segmentation = BAYSOR_RUN.out.segmentation
    ch_just_segmentation = ch_segmentation.map {
        _meta, segmentation -> return [ segmentation ]
    }
    ch_polygons2d = BAYSOR_RUN.out.polygons2d
    ch_htmls      = BAYSOR_RUN.out.htmls

    // run xeniumranger import-segmentation
    XENIUMRANGER_IMPORT_SEGMENTATION (
        ch_bundle_path,
        [],
        [],
        [],
        ch_just_segmentation,
        ch_polygons2d,
        "microns"
    )
    ch_versions = ch_versions.mix( XENIUMRANGER_IMPORT_SEGMENTATION.out.versions )

    ch_redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.bundle

    emit:

    segmentation     = ch_segmentation        // channel: [ val(meta), ["segmentation.csv"] ]
    polygons2d       = ch_polygons2d          // channel: [ ["segmentation_polygons_2d.json"] ]
    htmls            = ch_htmls               // channel: [ ["*.html"] ]

    redefined_bundle = ch_redefined_bundle    // channel: [ val(meta), "redefined-xenium-bundle" ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
