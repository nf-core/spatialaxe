//
// Run segger create_dataset, train and predict modules & parquet_to_csv
//

include { SEGGER2XR                        } from '../../../modules/local/utility/segger2xr/main'
include { SEGGER_TRAIN                     } from '../../../modules/local/segger/train/main'
include { SEGGER_PREDICT                   } from '../../../modules/local/segger/predict/main'
include { SEGGER_CREATE_DATASET            } from '../../../modules/local/segger/create_dataset/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow SEGGER_CREATE_TRAIN_PREDICT {

    take:

    ch_basedir             // channel: [ val(meta), [ "basedir" ] ]
    ch_transcripts_parquet // channel: [ val(meta), [bundle + "/transcripts.parquet"]]
    ch_bundle              // channel: [ val(meta), ["path-to-xenium-bundle"] ]

    main:

    ch_versions           = Channel.empty()
    ch_redefined_bundle   = Channel.empty()
    ch_segger_transcripts = Channel.empty()

    // create dataset
    SEGGER_CREATE_DATASET ( ch_basedir )
    ch_versions = ch_versions.mix ( SEGGER_CREATE_DATASET.out.versions )

    // train a model with the dataset created
    SEGGER_TRAIN ( SEGGER_CREATE_DATASET.out.datasetdir )
    ch_versions = ch_versions.mix ( SEGGER_TRAIN.out.versions )

    // run prediction with the trained models
    ch_just_trained_models = SEGGER_TRAIN.out.trained_models.map {
                _meta, models -> return [ models ]
    }
    ch_just_transcripts_parquet = ch_transcripts_parquet.map {
                _meta, transcripts -> return [ transcripts ]
    }
    SEGGER_PREDICT (
        SEGGER_CREATE_DATASET.out.datasetdir,
        ch_just_trained_models,
        ch_just_transcripts_parquet
    )
    ch_versions = ch_versions.mix ( SEGGER_PREDICT.out.versions )

    // convert parquet to XR compatible form
    SEGGER2XR ( SEGGER_PREDICT.out.transcripts )
    ch_versions = ch_versions.mix( SEGGER2XR.out.versions )

    ch_segger_transcripts = SEGGER2XR.out.transcripts_parquet.map {
        _meta, transcripts -> return [ transcripts ]
    }

    // replace transcripts.parquet in xenium bundle
    ch_updated_bundle = ch_bundle.map { fobj ->
        if (fobj.name == 'transcripts.parquet') {
            ch_segger_transcripts.val
        } else {
            fobj
        }
    }


    // run xeniumranger import-segmentation
    XENIUMRANGER_IMPORT_SEGMENTATION (
            ch_updated_bundle,
            [],
            [],
            [],
            [],
            [],
            "pixel"
        )
        ch_redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.bundle

        ch_versions = ch_versions.mix ( XENIUMRANGER_IMPORT_SEGMENTATION.out.versions )


    emit:

    datasetdir         = SEGGER_CREATE_DATASET.out.datasetdir // channel: [ val(meta), [ datasetdir ] ]
    trained_models     = SEGGER_TRAIN.out.trained_models      // channel: [ val(meta), [ trained_models ] ]
    benchmarks         = SEGGER_PREDICT.out.benchmarks        // channel: [ val(meta), [ benchmarks ] ]
    segger_transcripts = ch_segger_transcripts                // channel: [ [ transcripts.parquet ] ]

    redefined_bundle   = ch_redefined_bundle                  // channel: [ val(meta), ["redefined-xenium-bundle"] ]

    versions           = ch_versions                          // channel: [ versions.yml ]
}
