//
// Subworkflow with functionality specific to the nf-core/spatialxe pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args // array: List of positional nextflow CLI args
    outdir            // string: The output directory where the results will be saved
    input             // string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()
    log.info "INFO Input params validated  ✅ "

    //
    // Create channel from input file provided through params.input
    //
    try {

        Channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, bundle, image -> return [ [id: meta.id], bundle, image ]
        }
        .set { ch_samplesheet }

        log.info "INFO Samplesheet validated   ✅ "

    } catch (Exception e) {

        log.error "❌ Samplesheet validation failed: ${e.message}"
        exit 1
    }


    //
    // Check and validate xenium bundle
    //
    if ( !workflow.profile.contains('test')) {
        validateXeniumBundle(ch_samplesheet)
    }


    emit:

    samplesheet = ch_samplesheet
    versions    = ch_versions


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "❌ Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {

    // check if the segmentation method provided is valid for a mode
    if ( params.mode == 'image' && params.method ) {
        if ( !params.image_seg_methods.contains(params.method) ) {
            log.error "❌ Error: Invalid segmentation method: ${params.method} provided for the `image` based mode. Options: ${params.image_seg_methods}"
            exit 1
        }
    }

    if ( params.mode == 'coordinate' && params.method ) {
        if ( !params.transcript_seg_methods.contains(params.method) ) {
                log.error "❌ Error: Invalid segmentation method: `${params.method}` provided for the `coordinate` based mode. Options: ${params.transcript_seg_methods}"
                exit 1
        }
    }

    // check if --relabel_genes is true but --gene_panel is not provided
    if ( params.relabel_genes && !params.gene_panel ) {
        log.warn "⚠️  Relabel genes is enabled, but gene panel is not provided with the `--gene_panel`. Using `gene_panel.json` in the xenium bundle."
    }

    // check if --relabel_genes is true but --gene_panel is not provided
    if ( params.gene_panel && !params.relabel_genes ) {
        log.warn "⚠️  Gene panel provided, but relabel genes is disabled. Using `gene_panel.json` only to generate metadata."
    }

    // check if segmentation method is xeniumranger and nucleus_ony_segmentation is enabled
    if ( params.method == 'xeniumranger' && !params.nucleus_segmentation_only ) {
        log.warn "⚠️  Nucleus segmentation is disabled. Running xeniumranger resegment module to redefine xenium bundle without nucleus segmentation."
        log.warn "⚠️  Use --nucleus_segmentation_only to enable nucleus segmentation to redefine xenium bundle with import-segmentation module."
    }

    // check if segmentation mask is provided in image mode and baysor method
    if ( params.mode == 'image' && params.method == 'baysor' )
        if (!params.segmentation_mask ) {
        log.error "❌  Error: Missing segmentation mask with `--segmentation_mask` when pipeline is run in ${params.mode} and with the ${params.method}."
    }

}

//
// Check and validate xenium bundle
//
def validateXeniumBundle(ch_samplesheet) {

    // define xenium bundle directory structure
    def xenium_bundle = [
        "analysis.tar.gz",
        "analysis.zarr.zip",
        "analysis_summary.html",
        "cell_boundaries.csv.gz",
        "cell_boundaries.parquet",
        "cell_feature_matrix.h5",
        "cell_feature_matrix.tar.gz",
        "cell_feature_matrix.zarr.zip",
        "cells.csv.gz",
        "cells.parquet",
        "cells.zarr.zip",
        "experiment.xenium",
        "gene_panel.json",
        "metrics_summary.csv",
        "morphology.ome.tif",
        "morphology_focus/",
        "nucleus_boundaries.csv.gz",
        "nucleus_boundaries.parquet",
        "transcripts.parquet",
        "transcripts.zarr.zip"
    ]

    // get bundle path
    def ch_bundle_path = ch_samplesheet.map {
        _meta, bundle, _image ->
        def bundle_path = file (
            bundle.toString().replaceFirst(/\/$/, ''),
        )
        return bundle_path
    }

    // check if the path exists
    if ( !ch_bundle_path.map { it.exists() } ) {
        error "❌ Error: Xenium bundle path not found. Check if the path provided in the samplesheet exists."
        exit 1
    }

    // if the path exists, check for the presence of xenium files
    if ( ch_bundle_path.map { it.exists() } ) {

        ch_bundle_path.map { path ->
            def missing_files = []

            def allExist = xenium_bundle.every { filename ->
                def fullPath = file("${path}/${filename}")
                if (!fullPath.exists()) {
                    missing_files.add(filename)
                    return false
                }
                    return true
            }

            if (!allExist) {
                log.error "❌ Missing file(s) at bundle path provided in the samplesheet: ${missing_files}"
                exit 1
            }
        }
    }
    log.info "INFO Xenium bundle validated ✅ \n"
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

