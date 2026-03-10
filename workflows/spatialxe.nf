/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// multiqc
include { MULTIQC                                          } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_PRE_XR_RUN                    } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_POST_XR_RUN                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMultiqc                             } from '../subworkflows/nf-core/utils_nfcore_pipeline'

// nf-core functionality
include { softwareVersionsToYAML                           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                           } from '../subworkflows/local/utils_nfcore_spatialxe_pipeline'
include { paramsSummaryMap                                 } from 'plugin/nf-schema'

// nf-core modules
include { UNTAR                                            } from '../modules/nf-core/untar/main'

// coordinate-based segmentation subworklfows
include { SEGGER_CREATE_TRAIN_PREDICT                      } from '../subworkflows/local/segger_create_train_predict/main'
include { PROSEG_PRESET_PROSEG2BAYSOR                      } from '../subworkflows/local/proseg_preset_proseg2baysor/main'
include { PROSEG_PRESET_PROSEG2BAYSOR_TILED                } from '../subworkflows/local/proseg_preset_proseg2baysor_tiled/main'
include { BAYSOR_GENERATE_PREVIEW                          } from '../subworkflows/local/baysor_generate_preview/main'
include { BAYSOR_RUN_TRANSCRIPTS_PARQUET                   } from '../subworkflows/local/baysor_run_transcripts_parquet/main'

// image-based segmentation subworklfows
include { BAYSOR_RUN_PRIOR_SEGMENTATION_MASK               } from '../subworkflows/local/baysor_run_prior_segmentation_mask/main'
include { CELLPOSE_RESOLIFT_MORPHOLOGY_OME_TIF             } from '../subworkflows/local/cellpose_resolift_morphology_ome_tif/main'
include { CELLPOSE_BAYSOR_IMPORT_SEGMENTATION              } from '../subworkflows/local/cellpose_baysor_import_segmentation/main'
include { STARDIST_RESOLIFT_MORPHOLOGY_OME_TIF             } from '../subworkflows/local/stardist_resolift_morphology_ome_tif/main'
include { XENIUMRANGER_RESEGMENT_MORPHOLOGY_OME_TIF        } from '../subworkflows/local/xeniumranger_resegment_morphology_ome_tif/main'

// segmentation-free subworkflows
include { BAYSOR_GENERATE_SEGFREE                          } from '../subworkflows/local/baysor_generate_segfree/main'
include { FICTURE_PREPROCESS_MODEL                         } from '../subworkflows/local/ficture_preprocess_model/main'

// xeniumranger subworkflows
include { XENIUMRANGER_RELABEL_RESEGMENT                   } from '../subworkflows/local/xeniumranger_relabel_resegment/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION_REDEFINE_BUNDLE } from '../subworkflows/local/xeniumranger_import_segmentation_redefine_bundle/main'

// spatialdata subworkflows
include { SPATIALDATA_WRITE_META_MERGE                     } from '../subworkflows/local/spatialdata_write_meta_merge/main'

// TODO qc layer subworkflows
include { OPT_FLIP_TRACK_STAT                              } from '../subworkflows/local/opt_flip_track_stat/main'

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

    ch_versions = Channel.empty()

    ch_input = Channel.empty()
    ch_config = Channel.empty()
    ch_features = Channel.value([])
    ch_raw_bundle = Channel.empty()
    ch_gene_panel = Channel.empty()
    ch_qc_reports = Channel.empty()
    ch_bundle_path = Channel.empty()
    ch_preview_html = Channel.empty()
    ch_exp_metadata = Channel.empty()
    ch_gene_synonyms = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_multiqc_report = Channel.empty()
    ch_qupath_polygons = Channel.empty()
    ch_morphology_image = Channel.empty()
    ch_redefined_bundle = Channel.empty()
    ch_coordinate_space = Channel.empty()
    ch_panel_probes_fasta = Channel.empty()
    ch_transcripts_parquet = Channel.empty()
    ch_reference_annotations = Channel.empty()
    ch_multiqc_pre_xr_report = Channel.empty()
    ch_multiqc_post_xr_report = Channel.empty()


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - DATA STAGING
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // TODO: Replace with params.test_data_mode for robustness
    if (workflow.profile.contains('test')) {

        // get sample, xenium bundle and image path
        ch_input_untar = ch_samplesheet.map { meta, bundle, _image ->
            return [meta, bundle]
        }

        // get testdata
        UNTAR(ch_input_untar)

        ch_untar_outs = UNTAR.out.untar.map { meta, bundle ->
            return [meta, bundle.toString()]
        }

        ch_samplesheet
            .combine(ch_untar_outs, by: 0)
            .map { meta, _url, image, test_bundle ->
                return [meta, test_bundle, image]
            }
            .set { ch_input }
    }
    else {
        // for all other profile runs

        // check if samples are buffered
        if (params.buffer_samples) {
            ch_input = ch_samplesheet.buffer(size: params.buffer_size).map
            { buffered_sample ->
                def (meta, bundle, tif) = buffered_sample[0]
                tuple(meta, bundle, tif)
            }
        }
        else {
            ch_input = ch_samplesheet
        }
    }

    // path to bundle input
    ch_bundle_path = ch_input.map { meta, bundle, _image ->

        def bundle_path = file(bundle)
        if( !bundle_path.exists() ) {
            log.error("❌ Check if the path to the xenium bundle exists.")
            exit(1)
        }
        return [meta, bundle]
    }

    // get transcript.parquet from the xenium bundle
    ch_transcripts_parquet = ch_input.map { meta, bundle, _image ->
        def transcripts_parquet = file(
            bundle.toString().replaceFirst(/\/$/, '') + "/transcripts.parquet",
            checkIfExists: true
        )
        return [meta, transcripts_parquet]
    }

    // get morphology focus image from the xenium bundle (single 2D plane)
    // supports all Xenium versions:
    //   v2/v3: morphology_focus/morphology_focus_0000.ome.tif
    //   v4+:   morphology_focus/ch0000_dapi.ome.tif
    //   v1.x:  morphology_focus.ome.tif (single file at bundle root)
    //   fallback: morphology.ome.tif (multi-Z stack, not ideal for Cellpose)
    ch_morphology_image = ch_input.map { meta, bundle, image ->
        def morphology_img
        if (image) {
            morphology_img = file(image)
        } else {
            def bundle_path = bundle.toString().replaceFirst(/\/$/, '')
            def focus_v3 = file("${bundle_path}/morphology_focus/morphology_focus_0000.ome.tif")
            def focus_v4 = file("${bundle_path}/morphology_focus/ch0000_dapi.ome.tif")
            def focus_v1 = file("${bundle_path}/morphology_focus.ome.tif")
            if (focus_v3.exists()) {
                morphology_img = focus_v3
            } else if (focus_v4.exists()) {
                morphology_img = focus_v4
            } else if (focus_v1.exists()) {
                morphology_img = focus_v1
            } else {
                morphology_img = file("${bundle_path}/morphology.ome.tif", checkIfExists: true)
            }
        }
        return [meta, morphology_img]
    }

    // get experiment metdata - experiment.xenium
    ch_exp_metadata = ch_input.map { meta, bundle, _image ->
        def exp_metadata = file(
            bundle.toString().replaceFirst(/\/$/, '') + "/experiment.xenium",
            checkIfExists: true
        )
        return [meta, exp_metadata]
    }

    // get baysor xenium config
    ch_config = Channel.fromPath(
            "${projectDir}/assets/config/xenium.toml",
            checkIfExists: true
        )
        .flatten()

    // get segmentation mask if provided with --segmentation_mask for the baysor method
    if (params.segmentation_mask) {
        ch_segmentation_mask = Channel.fromPath(
                params.segmentation_mask,
                checkIfExists: true
            )
            .flatten()
    }

    // get a list of features if provided with the --features for the ficture method
    ch_features = params.features
        ? Channel.fromPath(params.features, checkIfExists: true).flatten()
        : Channel.value([])

    // get custom cellpose model if provided with the --cellpose_model for the cellpose method
    if (params.cellpose_model) {
        ch_cellpose_model = Channel.fromPath(
                params.cellpose_model,
                checkIfExists: true
            )
            .flatten()
    }

    // get panel probes fasta for off-target-probe tracking
    if (params.probes_fasta) {
        ch_panel_probes_fasta = Channel.fromPath(
                params.probes_fasta,
                checkIfExists: true
            )
            .flatten()
    }

    // get reference annotation files (gff,fa) for off-target-probe tracking
    if (params.reference_annotations) {
        ch_reference_annotations = Channel.fromPath(
                "${params.reference_annotations}/*.{fa,gff}".toString(),
                checkIfExists: true
            )
            .flatten()
    }

    // get gene synonyms for off-target-probe tracking
    if (params.gene_synonyms) {
        ch_gene_synonyms = Channel.fromPath(
                params.gene_synonyms,
                checkIfExists: true
            )
            .flatten()
    }

    // get qupath ploygons
    if (params.qupath_polygons) {
        ch_qupath_polygons = Channel.fromPath(
                "${params.qupath_polygons}/*.geojson",
                checkIfExists: true
            )
            .flatten()
    }

    // get gene_panel.json if provided with --gene_panel, sets relabel_genes to true
    def do_relabel = params.gene_panel ? true : params.relabel_genes
    if (params.gene_panel) {

        def gene_panel_file = file(params.gene_panel, checkIfExists: true)
        ch_gene_panel = ch_input.map { meta, _bundle, _image ->
            return [meta, gene_panel_file]
        }
    }
    else {

        // gene panel to use if only --relabel_genes is provided
        ch_gene_panel = ch_input.map { meta, bundle, _image ->
            def gene_panel = file(
                bundle.toString().replaceFirst(/\/$/, '') + "/gene_panel.json",
                checkIfExists: true
            )
            return [meta, gene_panel]
        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - RELABEL GENES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // run xr relabel if relabel_genes is true, check if gene_panel.json is provided
    if (do_relabel) {

        XENIUMRANGER_RELABEL_RESEGMENT(
            ch_bundle_path,
            ch_gene_panel,
        )
        ch_raw_bundle = XENIUMRANGER_RELABEL_RESEGMENT.out.redefined_bundle
    }
    else {
        ch_raw_bundle = ch_bundle_path
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - DATA PREVIEW
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // run baysor preview if `generate_preview ` is true
    if (params.mode == 'preview') {

        BAYSOR_GENERATE_PREVIEW(
            ch_transcripts_parquet,
            ch_config,
        )
        ch_preview_html = BAYSOR_GENERATE_PREVIEW.out.preview_html
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - XENIUMRANGER LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // run only xeniumranger import segmentation with changes xr specific params
    if (params.mode == 'image' && params.xeniumranger_only) {

        XENIUMRANGER_IMPORT_SEGMENTATION_REDEFINE_BUNDLE(
            ch_bundle_path
        )
        ch_redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION_REDEFINE_BUNDLE.out.redefined_bundle
        ch_coordinate_space = XENIUMRANGER_IMPORT_SEGMENTATION_REDEFINE_BUNDLE.out.coordinate_space
    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - IMAGE-BASED SEGMENTATION LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.mode == 'image') {

        // trigger the default image-based workflow if no method is specified
        if (!params.method) {

            CELLPOSE_BAYSOR_IMPORT_SEGMENTATION(
                ch_morphology_image,
                ch_bundle_path,
                ch_transcripts_parquet,
                ch_exp_metadata,
                ch_config,
            )
            ch_redefined_bundle = CELLPOSE_BAYSOR_IMPORT_SEGMENTATION.out.redefined_bundle
            ch_coordinate_space = CELLPOSE_BAYSOR_IMPORT_SEGMENTATION.out.coordinate_space
        }

        // run xeniumranger resegment with morphology_ome.tif
        if (params.method == 'xeniumranger') {

            XENIUMRANGER_RESEGMENT_MORPHOLOGY_OME_TIF(
                ch_bundle_path
            )
            ch_redefined_bundle = XENIUMRANGER_RESEGMENT_MORPHOLOGY_OME_TIF.out.redefined_bundle
            ch_coordinate_space = XENIUMRANGER_RESEGMENT_MORPHOLOGY_OME_TIF.out.coordinate_space
        }

        // run baysor run with morphology_ome.tif
        if (params.method == 'baysor') {

            if (params.segmentation_mask) {
                BAYSOR_RUN_PRIOR_SEGMENTATION_MASK(
                    ch_bundle_path,
                    ch_transcripts_parquet,
                    ch_segmentation_mask,
                    ch_config,
                )
            }
            ch_redefined_bundle = BAYSOR_RUN_PRIOR_SEGMENTATION_MASK.out.redefined_bundle
            ch_coordinate_space = BAYSOR_RUN_PRIOR_SEGMENTATION_MASK.out.coordinate_space
        }

        // run cellpose on the morphology_ome.tif
        if (params.method == 'cellpose') {

            CELLPOSE_RESOLIFT_MORPHOLOGY_OME_TIF(
                ch_morphology_image,
                ch_bundle_path,
            )
            ch_redefined_bundle = CELLPOSE_RESOLIFT_MORPHOLOGY_OME_TIF.out.redefined_bundle
            ch_coordinate_space = CELLPOSE_RESOLIFT_MORPHOLOGY_OME_TIF.out.coordinate_space
        }

        // run stardist on the morphology_ome.tif
        if (params.method == 'stardist') {

            STARDIST_RESOLIFT_MORPHOLOGY_OME_TIF(
                ch_morphology_image,
                ch_bundle_path,
            )
            ch_redefined_bundle = STARDIST_RESOLIFT_MORPHOLOGY_OME_TIF.out.redefined_bundle
            ch_coordinate_space = STARDIST_RESOLIFT_MORPHOLOGY_OME_TIF.out.coordinate_space
        }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - TRANSCRIPT-BASED SEGMENTATION LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.mode == 'coordinate') {

        // run proseg with transcripts.parquet if method = proseg or is not provided (default workflow)
        if (!params.method || params.method == 'proseg') {

            if (params.tiling) {
                PROSEG_PRESET_PROSEG2BAYSOR_TILED(
                    ch_bundle_path,
                    ch_transcripts_parquet,
                )
                ch_redefined_bundle = PROSEG_PRESET_PROSEG2BAYSOR_TILED.out.redefined_bundle
                ch_coordinate_space = PROSEG_PRESET_PROSEG2BAYSOR_TILED.out.coordinate_space
            } else {
                PROSEG_PRESET_PROSEG2BAYSOR(
                    ch_bundle_path,
                    ch_transcripts_parquet,
                )
                ch_redefined_bundle = PROSEG_PRESET_PROSEG2BAYSOR.out.redefined_bundle
                ch_coordinate_space = PROSEG_PRESET_PROSEG2BAYSOR.out.coordinate_space
            }
        }

        // run segger with transcripts.parquet
        if (params.method == 'segger') {

            SEGGER_CREATE_TRAIN_PREDICT(
                ch_bundle_path,
                ch_transcripts_parquet,
            )
            ch_redefined_bundle = SEGGER_CREATE_TRAIN_PREDICT.out.redefined_bundle
            ch_coordinate_space = SEGGER_CREATE_TRAIN_PREDICT.out.coordinate_space
        }

        // run baysor with transcripts.parquet (unified tiled/non-tiled subworkflow)
        if (params.method == 'baysor') {

            // Image-based prior (cellpose mask) requires non-tiled Baysor
            if ( params.baysor_tiling && params.baysor_prior == 'cellpose' ) {
                error "ERROR: baysor_prior='cellpose' (image-based) requires baysor_tiling=false. " +
                      "For tiled Baysor, use baysor_prior='cells' (column-based)."
            }

            ch_prior_mask = Channel.empty()

            BAYSOR_RUN_TRANSCRIPTS_PARQUET(
                ch_bundle_path,
                ch_transcripts_parquet,
                ch_morphology_image,
                ch_config,
                ch_prior_mask,
            )
            ch_redefined_bundle = BAYSOR_RUN_TRANSCRIPTS_PARQUET.out.redefined_bundle
            ch_coordinate_space = BAYSOR_RUN_TRANSCRIPTS_PARQUET.out.coordinate_space
        }
    }



    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - SPATIALDATA / METADATA LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // run spatialdata modules to generate sd objects in image or coordinate mode
    if (params.mode == 'image' || params.mode == 'coordinate') {

        SPATIALDATA_WRITE_META_MERGE(
            ch_bundle_path,
            ch_redefined_bundle,
            ch_coordinate_space,
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - QC LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // check to run the qc layer
    if (params.mode == 'qc' || params.run_qc) {

        if (params.offtarget_probe_tracking) {

            // run off-target probe tracking
            OPT_FLIP_TRACK_STAT(
                ch_panel_probes_fasta,
                ch_reference_annotations,
                ch_gene_synonyms,
            )
        }
    }


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - SEGMENTATION-FREE LAYER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.mode == 'segfree') {

        // trigger the default segfree workflow if no method or if the method is baysor
        if (!params.method || params.method == 'baysor') {

            BAYSOR_GENERATE_SEGFREE(
                ch_transcripts_parquet,
                ch_config,
            )
        }

        // run ficture with transcripts.parquet
        if (params.method == 'ficture') {

            FICTURE_PREPROCESS_MODEL(
                ch_transcripts_parquet,
                ch_features,
            )
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
            name: 'nf_core_' + 'spatialxe_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SPATIALXE - MultiQC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )

    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()

    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )

    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )

    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)

    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    if (params.mode == 'image' || params.mode == 'coordinate') {

        // get path to the raw bundle
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_bundle_path.map { _meta, bundle -> file(bundle) }.collect().ifEmpty([])
        )

        MULTIQC_PRE_XR_RUN (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            [],
        )
        ch_multiqc_pre_xr_report = MULTIQC_PRE_XR_RUN.out.report.toList()

        // get path to the redefined bundle
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_redefined_bundle.map { _meta, bundle -> file(bundle) }.collect().ifEmpty([])
        )

        MULTIQC_POST_XR_RUN (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            [],
        )
        ch_multiqc_post_xr_report = MULTIQC_POST_XR_RUN.out.report.toList()

    } else {

        // get path to the raw bundle
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_bundle_path.map { _meta, bundle -> file(bundle) }.collect().ifEmpty([])
        )


        // get the qc htmls if qc mode is run
        if (params.mode == 'qc' || params.run_qc) {

            // TODO collect all qc outs in a channel to be passed to multiqc
            ch_multiqc_files = ch_multiqc_files.mix(
                ch_qc_reports.map { _meta, qc_reports -> qc_reports }.collect().ifEmpty([])
            )

        }


        // get the preview html if preview mode is run
        if (params.mode == 'preview') {

            ch_multiqc_files = ch_multiqc_files.mix(
                ch_preview_html.map { _meta, preview_html -> preview_html }.collect().ifEmpty([])
            )

        }


        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            [],
        )
        ch_multiqc_report = MULTIQC.out.report.toList()

    }

    emit:
    multiqc_pre_xr_report  = ch_multiqc_pre_xr_report  // channel: /path/to/multiqc_report.html
    multiqc_post_xr_report = ch_multiqc_post_xr_report // channel: /path/to/multiqc_report.html
    multiqc_report         = ch_multiqc_report         // channel: /path/to/multiqc_report.html
    versions               = ch_versions               // channel: [ path(versions.yml) ]
}
