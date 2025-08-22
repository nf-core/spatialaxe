//
// Run ficture preprocess and model modules
//

include { FICTURE_PREPROCESS } from '../../../modules/local/ficture/preprocess/main'
include { FICTURE            } from '../../../modules/local/ficture/model/main'
include { PARQUET_TO_CSV     } from '../../../modules/local/spatialconverter/parquet_to_csv/main'


workflow FICTURE_PREPROCESS_MODEL {

    take:

    ch_transcripts_parquet // channel: [ val(meta), [ "transcripts.parquet" ] ]
    ch_features            // channel: [ ["features"] ]

    main:

    ch_versions = Channel.empty()

    // convert parquet to csv
    PARQUET_TO_CSV ( ch_transcripts_parquet, ".csv" )
    ch_versions = ch_versions.mix ( PARQUET_TO_CSV.out.versions )

    // run ficture preprocessing
    ch_transcripts = PARQUET_TO_CSV.out.transcripts_csv

    FICTURE_PREPROCESS ( ch_transcripts, ch_features )
    ch_versions = ch_versions.mix ( FICTURE_PREPROCESS.out.versions )

    // run the ficture wrapper pipeline
    ch_features_clean = Channel.empty()
    if ( params.features ) {
        ch_features_clean = FICTURE_PREPROCESS.out.features
    }
    FICTURE (
        FICTURE_PREPROCESS.out.transcripts,
        FICTURE_PREPROCESS.out.coordinate_minmax,
        ch_features_clean
    )
    ch_versions = ch_versions.mix( FICTURE.out.versions )

    emit:

    transcripts        = FICTURE_PREPROCESS.out.transcripts         // channel: [ val(meta), [ "*processed_transcripts.tsv.gz" ] ]
    coordinate_minmax  = FICTURE_PREPROCESS.out.coordinate_minmax   // channel: [ "*coordinate_minmax.tsv" ]
    features           = FICTURE_PREPROCESS.out.features            // channel: [ "*feature.clean.tsv.gz" ]

    results            = FICTURE.out.results                        // channel: [ val(meta), [ "results/** ] ]

    versions = ch_versions                                          // channel: [ versions.yml ]
}

