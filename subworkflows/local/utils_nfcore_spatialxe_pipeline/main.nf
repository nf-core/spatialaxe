//
// Subworkflow with functionality specific to the nf-core/spatialxe pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { completionEmail         } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary       } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification          } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

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
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/spatialxe ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { "    https://doi.org/${it.trim().replace('https://doi.org/', '')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/spatialxe/blob/master/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --mode <MODE> --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()
    log.info("✅ Pipeline parameters validated.")

    //
    // Create channel from input file provided through params.input
    //
    try {

        Channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
            .map { meta, bundle, image ->
                return [[id: meta.id], bundle, image]
            }
            .set { ch_samplesheet }

        log.info("✅ Samplesheet validated.")
    }
    catch (Exception e) {

        log.error("❌ Samplesheet validation failed: ${e.message}")
        exit(1)
    }


    //
    // Check and validate xenium bundle
    //
    if (!workflow.profile.contains('test')) {
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
        log.error("❌ Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
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

    // check if conda profile is provided
    if (workflow.profile.contains('conda')) {
        log.error("❌ Error: `nf-core/spatialxe` does not support running the pipeline with profile: conda ")
        exit(1)
    }

    // check if the samplesheet provided with the test config is assets/samplesheet.csv
    if (workflow.profile.contains('test') && !"${params.input}".endsWith("assets/samplesheet.csv")) {
        log.error("❌ Error: Use the samplesheet at: ${projectDir}/assets/samplesheet.csv with `--input` when running the pipeline in test profile.")
        exit(1)
    }

    // check if the segmentation method provided is valid for a mode
    if (params.mode == 'image' && params.method) {
        if (!params.image_seg_methods.contains(params.method)) {
            log.error("❌ Error: Invalid segmentation method: ${params.method} provided for the `image` based mode. Options: ${params.image_seg_methods}")
            exit(1)
        }
    }

    if (params.mode == 'coordinate' && params.method) {
        if (!params.transcript_seg_methods.contains(params.method)) {
            log.error("❌ Error: Invalid segmentation method: `${params.method}` provided for the `coordinate` based mode. Options: ${params.transcript_seg_methods}")
            exit(1)
        }
    }

    // check if --relabel_genes is true but --gene_panel is not provided
    if (params.relabel_genes && !params.gene_panel) {
        log.warn("⚠️  Relabel genes is enabled, but gene panel is not provided with the `--gene_panel`. Using `gene_panel.json` in the xenium bundle.")
    }

    // check if --relabel_genes is true but --gene_panel is not provided
    if (params.gene_panel && !params.relabel_genes) {
        log.warn("⚠️  Gene panel provided, but relabel genes is disabled. Using `gene_panel.json` only to generate metadata.")
    }

    // check if segmentation method is xeniumranger and nucleus_ony_segmentation is enabled
    if (params.method == 'xeniumranger' && !params.nucleus_segmentation_only) {
        log.warn("⚠️  Nucleus segmentation is disabled. Running xeniumranger resegment module to redefine xenium bundle without nucleus segmentation.")
        log.warn("⚠️  Use --nucleus_segmentation_only to enable nucleus segmentation to redefine xenium bundle with import-segmentation module.")
    }

    // check if segmentation mask is provided in image mode and baysor method
    if (params.mode == 'image' && params.method == 'baysor') {
        if (!params.segmentation_mask) {
            log.warn("⚠️  Missing segmentation mask with `--segmentation_mask` when pipeline is run in ${params.mode} and with the ${params.method}. Running in coordinate mode.")
        }
    }

    // check if required arguments are provided for off-target probe tracking
    if (!params.mode && params.offtarget_probe_tracking) {
        if(!params.probes_fasta || !params.reference_annotations || !params.gene_synonyms) {
            log.error("❌ Error: Missing required param(s) for off-target-proebe detection.")
            exit(1)
        }
        log.error("❌ Error: Use --mode qc and --offtraget_probe_tracking to run off-target probe tracking.")
        exit(1)
    }
}

//
// Check and validate xenium bundle
//
def validateXeniumBundle(ch_samplesheet) {

    // define xenium bundle directory structure - required files
    def bundle_required_files = [
        "cell_boundaries.csv.gz",
        "cell_boundaries.parquet",
        "cell_feature_matrix.h5",
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
        "transcripts.zarr.zip",
    ]

    // bundle optional files
    def bundle_optional_files = [
        "analysis.tar.gz",
        "analysis.zarr.zip",
        "analysis_summary.html"
    ]

    // get bundle path (keep raw string for remote-path detection)
    def ch_bundle_info = ch_samplesheet.map { _meta, bundle, _image ->
        def rawPath = bundle.toString().replaceFirst(/\/$/, '')
        def bundle_path = file(rawPath)
        return [rawPath, bundle_path]
    }

    // Skip file-level validation for remote paths (S3, GS, AZ) because
    // file().exists() is unreliable on cloud storage during initialization
    // (Fusion mounts s3://bucket as /bucket, breaking startsWith checks).
    // Files will be validated at task staging time instead.
    ch_bundle_info.map { rawPath, path ->
        if (rawPath.startsWith('s3://') || rawPath.startsWith('gs://') || rawPath.startsWith('az://')) {
            log.info("Skipping bundle file validation for remote path: ${rawPath}")
            return
        }

        def missing_required_files = []
        def missing_optional_files = []

        def requiredExist = bundle_required_files.every { filename ->
            def fullPath = file("${path}/${filename}")
            if (!fullPath.exists()) {
                missing_required_files.add(filename)
                return false
            }
            return true
        }
        if (!requiredExist) {
            log.error("❌ Missing file(s) at bundle path provided in the samplesheet: ${missing_required_files}")
            exit(1)
        }

        def optionalExist = bundle_optional_files.every { filename ->
            def fullPath = file("${path}/${filename}")
            if (!fullPath.exists()) {
                missing_optional_files.add(filename)
                return false
            }
            return true
        }
        if (!optionalExist) {
            log.warn("⚠️ Missing optional file(s) at bundle path provided in the samplesheet: ${missing_optional_files}")
        }

        log.info("✅ Xenium bundle validated.\n")
    }
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
        ".",
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
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
