//
// Run baysor create_dataset & preview
//

include { BAYSOR_PREVIEW        } from '../../../modules/local/baysor/preview/main'
include { BAYSOR_CREATE_DATASET } from '../../../modules/local/baysor/create_dataset/main'
include { PARQUET_TO_CSV        } from '../../../modules/local/spatialconverter/parquet_to_csv/main'

workflow BAYSOR_GENERATE_PREVIEW {

    take:

    ch_transcripts_parquet // channel: [ val(meta), ["path-to-transcripts.parquet"] ]
    ch_config              // channel: ["path-to-xenium.toml"]

    main:

    ch_versions      = Channel.empty()
    ch_preview_html  = Channel.empty()


    // run parquet to csv
    PARQUET_TO_CSV ( ch_transcripts_parquet, ".csv" )
    ch_versions = ch_versions.mix ( PARQUET_TO_CSV.out.versions )

    // generate randomised sample data
    BAYSOR_CREATE_DATASET ( PARQUET_TO_CSV.out.transcripts_csv, 0.3 )
    ch_versions = ch_versions.mix ( BAYSOR_CREATE_DATASET.out.versions )

    // run baysor preview if param - generate_preview is true
    BAYSOR_PREVIEW (
        BAYSOR_CREATE_DATASET.out.sampled_transcripts,
        ch_config
    )
    ch_versions = ch_versions.mix ( BAYSOR_PREVIEW.out.versions )

    ch_preview_html = BAYSOR_PREVIEW.out.preview_html

    emit:

    preview_html     = ch_preview_html   // channel: [ val(meta), ["preview.html"] ]

    versions         = ch_versions       // channel: [ versions.yml ]
}
