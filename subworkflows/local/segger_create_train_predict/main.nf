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

    ch_bundle              // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    ch_transcripts_parquet // channel: [ val(meta), [bundle + "/transcripts.parquet"]]


    main:

    ch_versions           = Channel.empty()

    ch_updated_bundle     = Channel.empty()
    ch_redefined_bundle   = Channel.empty()
    ch_segger_transcripts = Channel.empty()
    ch_coordinate_space   = Channel.value("pixels")

    // create dataset
    SEGGER_CREATE_DATASET ( ch_bundle )
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

    ch_segger_transcripts_parquet = SEGGER2XR.out.transcripts_parquet.map {
        _meta, transcripts -> return [ transcripts ]
    }


    // swap transscripts.parquet with segger transcripts
    ch_updated_bundle = ch_bundle.map { _meta, fileobj ->
        fileobj == "transcripts.parquet" ? ch_segger_transcripts_parquet : fileobj
    }


    // run xeniumranger import-segmentation
    cells = ch_updated_bundle.map { _meta, bundle ->
        return [ bundle + "/cells.zarr.zip" ]
    }

    if ( params.nucleus_segmentation_only ) {

        XENIUMRANGER_IMPORT_SEGMENTATION (
            ch_updated_bundle,
            [],
            cells,
            cells,
            [],
            [],
            ch_coordinate_space
        )
        ch_redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.bundle

        ch_versions = ch_versions.mix ( XENIUMRANGER_IMPORT_SEGMENTATION.out.versions )

    } else {

        XENIUMRANGER_IMPORT_SEGMENTATION (
            ch_updated_bundle,
            [],
            [],
            cells,
            [],
            [],
            ch_coordinate_space
        )
        ch_redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.bundle

        ch_versions = ch_versions.mix ( XENIUMRANGER_IMPORT_SEGMENTATION.out.versions )

    }

    emit:

    datasetdir         = SEGGER_CREATE_DATASET.out.datasetdir // channel: [ val(meta), [ datasetdir ] ]
    trained_models     = SEGGER_TRAIN.out.trained_models      // channel: [ val(meta), [ trained_models ] ]
    benchmarks         = SEGGER_PREDICT.out.benchmarks        // channel: [ val(meta), [ benchmarks ] ]
    segger_transcripts = ch_segger_transcripts                // channel: [ [ transcripts.parquet ] ]

    coordinate_space   = ch_coordinate_space                  // channel: [ ["pixels"] ]

    redefined_bundle   = ch_redefined_bundle                  // channel: [ val(meta), ["redefined-xenium-bundle"] ]

    versions           = ch_versions                          // channel: [ versions.yml ]
}
