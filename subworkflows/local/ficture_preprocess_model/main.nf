//
// Run ficture preprocess and model modules
//

include { FICTURE_PREPROCESS } from '../../../modules/local/ficture/preprocess/main'
include { FICTURE            } from '../../../modules/local/ficture/model/main'
include { PARQUET_TO_CSV     } from '../../../modules/local/utility/parquet_to_csv/main'



workflow FICTURE_PREPROCESS_MODEL {
    take:
    ch_transcripts_file // channel: [ val(meta), [ "transcripts.parquet" ] ]
    ch_features            // channel: [ ["features"] ]
    features               // value: path to features list (or null)

    main:

    // convert parquet to csv
    PARQUET_TO_CSV(ch_transcripts_file, ".csv")

    // run ficture preprocessing
    ch_transcripts = PARQUET_TO_CSV.out.transcripts_csv

    FICTURE_PREPROCESS(ch_transcripts, ch_features)

    // run the ficture wrapper pipeline
    ch_features_clean = features ? FICTURE_PREPROCESS.out.features : channel.value([])
    FICTURE(
        FICTURE_PREPROCESS.out.transcripts,
        FICTURE_PREPROCESS.out.coordinate_minmax,
        ch_features_clean,
    )
    emit:
    transcripts       = FICTURE_PREPROCESS.out.transcripts // channel: [ val(meta), [ "*processed_transcripts.tsv.gz" ] ]
    coordinate_minmax = FICTURE_PREPROCESS.out.coordinate_minmax // channel: [ "*coordinate_minmax.tsv" ]
    features          = FICTURE_PREPROCESS.out.features // channel: [ "*feature.clean.tsv.gz" ]
    results           = FICTURE.out.results // channel: [ val(meta), [ "results/** ] ]
}
