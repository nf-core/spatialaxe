/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// multiqc
include { MULTIQC              } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'

// nf-core functionality
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_spatialxe_pipeline'
include { paramsSummaryMap       } from 'plugin/nf-schema'

// nf-core modules
include { UNTAR } from '../modules/nf-core/untar/main'

// coordinate-based segmentation subworklfows
include { SEGGER_CREATE_TRAIN_PREDICT    } from '../subworkflows/local/segger_create_train_predict/main'
include { PROSEG_PRESET_PROSEG2BAYSOR    } from '../subworkflows/local/proseg_preset_proseg2baysor/main'
include { BAYSOR_GENERATE_PREVIEW        } from '../subworkflows/local/baysor_generate_preview/main'
include { BAYSOR_RUN_TRANSCRIPTS_PARQUET } from '../subworkflows/local/baysor_run_transcripts_parquet/main'

// image-based segmentation subworklfows
include { BAYSOR_RUN_PRIOR_SEGMENTATION_MASK        } from '../subworkflows/local/baysor_run_prior_segmentation_mask/main'
include { CELLPOSE_RESOLIFT_MORPHOLOGY_OME_TIF      } from '../subworkflows/local/cellpose_resolift_morphology_ome_tif/main'
include { CELLPOSE_BAYSOR_IMPORT_SEGMENTATION       } from '../subworkflows/local/cellpose_baysor_import_segmentation/main'
include { XENIUMRANGER_RESEGMENT_MORPHOLOGY_OME_TIF } from '../subworkflows/local/xeniumranger_resegment_morphology_ome_tif/main'

// segmentation-free subworkflows
include { BAYSOR_GENERATE_SEGFREE  } from '../subworkflows/local/baysor_generate_segfree/main'
include { FICTURE_PREPROCESS_MODEL } from '../subworkflows/local/ficture_preprocess_model/main'

// xeniumranger subworkflows
include { XENIUMRANGER_RELABEL_RESEGMENT                   } from '../subworkflows/local/xeniumranger_relabel_resegment/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION_REDEFINE_BUNDLE } from '../subworkflows/local/xeniumranger_import_segmentation_redefine_bundle/main'

// spatialdata subworkflows
include { SPATIALDATA_WRITE_META_MERGE } from '../subworkflows/local/spatialdata_write_meta_merge/main'

// TODO qc layer subworkflows


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPATIALXE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - GENERATE INPUT CHANNELS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_versions            = Channel.empty()

    ch_input               = Channel.empty()
    ch_bundle              = Channel.empty()
    ch_config              = Channel.empty()
    ch_features            = Channel.empty()
    ch_raw_bundle          = Channel.empty()
    ch_gene_panel          = Channel.empty()
    ch_bundle_path         = Channel.empty()
    ch_multiqc_files       = Channel.empty()
    ch_morphology_image    = Channel.empty()
    ch_redefined_bundle    = Channel.empty()
    ch_transcripts_parquet = Channel.empty()


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - DATA STAGING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if ( workflow.profile.contains('test') ) {

        // get sample, xenium bundle and image path
        ch_input_untar = ch_samplesheet.map { meta, bundle, _image ->
            return [ meta, bundle ]
        }

        // get testdata
        UNTAR ( ch_input_untar )
        ch_versions = ch_versions.mix ( UNTAR.out.versions )

        ch_untar_outs = UNTAR.out.untar.map { _meta, bundle ->
            return [ bundle.toString() ]
        }

        ch_samplesheet.combine(ch_untar_outs)
            .map { meta, _url, image, test_bundle ->
                return [meta, test_bundle, image]
            }
            .set { ch_input }

    } else {

        // for all other profile runs
        ch_input = ch_samplesheet
    }

    // path to bundle input
    ch_bundle_path = ch_input.map { meta, bundle, _image ->
        return [ meta, bundle ]
    }

    // get transcript.parquet from the xenium bundle
    ch_transcripts_parquet = ch_input.map { meta, bundle, _image ->
        def transcripts_parquet = file (
            bundle.toString().replaceFirst(/\/$/, '') + "/transcripts.parquet",
            checkIfExists: true
        )
        return [ meta, transcripts_parquet ]
    }

    // get morphology.ome.tif from the xenium bundle
    ch_morphology_image = ch_input.map { meta, bundle, image ->
        def morphology_img = image ? file(image) : file (
            bundle.toString().replaceFirst(/\/$/, '') + "/morphology.ome.tif",
            checkIfExists: true
        )
        return [ meta, morphology_img ]
    }

    // get baysor xenium config
    ch_config = Channel.fromPath (
        "${projectDir}/assets/config/xenium.toml",
        checkIfExists: true
    )

    // get segmentation mask if provided with --segmentation_mask for the baysor method
    if ( params.segmentation_mask ) {
        ch_segmentation_mask = Channel.fromPath (
            params.segmentation_mask, checkIfExists: true
        )
    }

    // get a list of features if provided with the --features for the ficture method
    if ( params.features ) {
        ch_features = Channel.fromPath (
            params.features,
            checkIfExists: true
        )
    }

    // get custom cellpose model if provided with the --cellpose_model for the cellpose method
    if ( params.cellpose_model ) {
        ch_features = Channel.fromPath (
            params.cellpose_model,
            checkIfExists: true
        )
    }

    // get gene_panel.json if provided with --gene_panel, sets relabel_genes to true
    if (( params.gene_panel )) {

        params.relabel_genes = true
        ch_gene_panel = Channel.fromPath (
            params.gene_panel,
            checkIfExists: true
        )

    } else {

        // gene panel to use if only --relabel_genes is provided
        ch_gene_panel = ch_input.map { meta, bundle, _image ->
            def gene_panel = file (
                bundle.toString().replaceFirst(/\/$/, '') + "/gene_panel.json",
                checkIfExists: true
            )
            return [ meta, gene_panel ]
        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - RELABEL GENES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // run xr relabel if relabel_genes is true, check if gene_panel.json is provided
    if ( params.relabel_genes ) {

        XENIUMRANGER_RELABEL_RESEGMENT (
            ch_bundle_path,
            ch_gene_panel
        )
        ch_raw_bundle = XENIUMRANGER_RELABEL_RESEGMENT.out.redefined_bundle

    } else {
        ch_raw_bundle = ch_bundle
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - DATA PREVIEW
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // run baysor preview if `generate_preview ` is true
    if ( params.mode == 'preview' ) {

        BAYSOR_GENERATE_PREVIEW (
            ch_transcripts_parquet,
            ch_config
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - IMAGE-BASED SEGMENTATION LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if ( params.mode == 'image' ) {

        // trigger the default image-based workflow if no method is specified
        if ( !params.method ) {

            CELLPOSE_BAYSOR_IMPORT_SEGMENTATION (
                ch_morphology_image,
                ch_bundle_path,
                ch_transcripts_parquet,
                ch_config
            )
            ch_redefined_bundle = CELLPOSE_BAYSOR_IMPORT_SEGMENTATION.out.redefined_bundle
        }

        // check it the provided method is part of the methods list
        if ( params.method in params.image_seg_methods ) {

            // run xeniumranger resegment with morphology_ome.tif
            if ( params.method == 'xeniumranger' ) {

                XENIUMRANGER_RESEGMENT_MORPHOLOGY_OME_TIF (
                    ch_bundle_path
                )
                ch_redefined_bundle = XENIUMRANGER_RESEGMENT_MORPHOLOGY_OME_TIF.out.redefined_bundle
            }

            // run baysor run with morphology_ome.tif
            if ( params.method == 'baysor' ) {

                if ( params.segmentation_mask ) {
                    BAYSOR_RUN_PRIOR_SEGMENTATION_MASK (
                        ch_bundle_path,
                        ch_transcripts_parquet,
                        ch_segmentation_mask,
                        ch_config
                    )
                } else {
                    BAYSOR_RUN_PRIOR_SEGMENTATION_MASK (
                        ch_bundle_path,
                        ch_transcripts_parquet,
                        [],
                        ch_config
                    )
                }
                ch_redefined_bundle = BAYSOR_RUN_PRIOR_SEGMENTATION_MASK.out.redefined_bundle
            }

            // run cellpose on the morphology_ome.tif
            if ( params.method == 'cellpose' ) {

                CELLPOSE_RESOLIFT_MORPHOLOGY_OME_TIF (
                    ch_morphology_image,
                    ch_bundle_path
                )
                ch_redefined_bundle = CELLPOSE_RESOLIFT_MORPHOLOGY_OME_TIF.out.redefined_bundle
            }

        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - TRANSCRIPT-BASED SEGMENTATION LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if ( params.mode == 'coordinate' ) {

        // trigger the default transcripts-based workflow if no method is specified
        if ( !params.method ) {

            PROSEG_PRESET_PROSEG2BAYSOR (
                ch_bundle_path,
                ch_transcripts_parquet
            )
            ch_redefined_bundle = PROSEG_PRESET_PROSEG2BAYSOR.out.redefined_bundle

        }

        // check it the provided method is part of the methods list
        if ( params.method in params.transcript_seg_methods ) {

            // run proseg with transcripts.parquet
            if ( params.method == 'proseg') {

                PROSEG_PRESET_PROSEG2BAYSOR (
                    ch_bundle_path,
                    ch_transcripts_parquet
                )
                ch_redefined_bundle = PROSEG_PRESET_PROSEG2BAYSOR.out.redefined_bundle

            }

            // run segger with transcripts.parquet
            if ( params.method == 'segger' ) {

                SEGGER_CREATE_TRAIN_PREDICT (
                    ch_bundle_path,
                    ch_transcripts_parquet
                )

            }

            // run baysor with transcripts.parquet
            if ( params.method == 'baysor' ) {

                BAYSOR_RUN_TRANSCRIPTS_PARQUET (
                    ch_bundle_path,
                    ch_transcripts_parquet,
                    ch_config
                )
                ch_redefined_bundle = BAYSOR_RUN_TRANSCRIPTS_PARQUET.out.redefined_bundle
            }

        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - XENIUMRANGER LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // run only xeniumranger import segmentation with changes xr specific params
    if  ( params.xeniumranger_only ) {

        XENIUMRANGER_IMPORT_SEGMENTATION_REDEFINE_BUNDLE (
            ch_bundle_path
        )
        ch_redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION_REDEFINE_BUNDLE.out.redefined_bundle
    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - SPATIALDATA / METADATA LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // run spatialdata modules to generate sd objects in image or coordinate mode
    if ( params.mode == 'image' || params.mode == 'coordinate' ) {

        SPATIALDATA_WRITE_META_MERGE (
            ch_bundle_path,
            ch_redefined_bundle
        )

    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - QC LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - SEGMENTATION-FREE LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if ( params.mode == 'segfree' ) {

        // trigger the default segfree workflow if no method is specified
        if ( !params.method ) {

            BAYSOR_GENERATE_SEGFREE (
                ch_transcripts_parquet,
                ch_config
            )
        }

        // check it the provided method is part of the methods list
        if ( params.method in params.segfree_methods ) {

            // run baysor with transcripts.parquet
            if ( params.method == 'baysor' ) {

                BAYSOR_GENERATE_SEGFREE (
                    ch_transcripts_parquet,
                    ch_config
                )
            }

            // run ficture with transcripts.parquet
            if ( params.method == 'ficture' ) {

                FICTURE_PREPROCESS_MODEL (
                    ch_transcripts_parquet,
                    ch_features
                )
            }
        }
    }



    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - COLLATE & SAVE SOFTWARE VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'spatialxe_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - MultiQC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_multiqc_config        = Channel.fromPath (
        "$projectDir/assets/multiqc_config.yml",
        checkIfExists: true
    )

    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath ( params.multiqc_config, checkIfExists: true ) :
        Channel.empty()

    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath( params.multiqc_logo, checkIfExists: true ) :
        Channel.empty()

    summary_params      = paramsSummaryMap (
        workflow, parameters_schema: "nextflow_schema.json"
    )

    ch_workflow_summary = Channel.value( paramsSummaryMultiqc( summary_params ) )

    ch_multiqc_files = ch_multiqc_files.mix (
        ch_workflow_summary.collectFile( name: 'workflow_summary_mqc.yaml' )
    )

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file( params.multiqc_methods_description, checkIfExists: true ) :
        file( "$projectDir/assets/methods_description_template.yml", checkIfExists: true )

    ch_methods_description                = Channel.value (
        methodsDescriptionText ( ch_multiqc_custom_methods_description )
    )

    ch_multiqc_files = ch_multiqc_files.mix ( ch_collated_versions )

    ch_multiqc_files = ch_multiqc_files.mix (
        ch_methods_description.collectFile (
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:

    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
