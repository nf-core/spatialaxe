//
// IMAGE_QC: image-based quality control for a Xenium bundle.
//
// Runs the image QC analysis (focus / SNR / morphology metrics + figures) and
// renders an HTML report from the analysis outputs with the nf-core
// QUARTONOTEBOOK module. Versions are reported via the `versions` topic channel
// by each module (collected centrally by the pipeline), so this subworkflow
// does not thread versions through emit.
//

include { IMAGE_QC_ANALYSIS as ANALYSIS } from '../../../modules/local/image_qc/main'
include { QUARTONOTEBOOK as REPORT      } from '../../../modules/nf-core/quartonotebook/main'

workflow IMAGE_QC {
    take:
    ch_input // channel: [ val(meta), val(parameters), path(input_files) ]
    notebook // path: the image QC report .qmd

    main:
    // ROI threshold config: user-provided path, else the bundled default.
    ch_roi_yaml = channel.fromPath(
        params.roi_image_qc_thresholds_yaml ?: "${projectDir}/conf/roi_image_qc_thresholds.yaml",
        checkIfExists: true,
    )

    // Single process handles both ROI-based and cell-level QC
    ANALYSIS(ch_input, ch_roi_yaml.first())

    // Render the report with QUARTONOTEBOOK. Its four inputs are separate
    // channels paired by emission order, so all per-sample channels are derived
    // from the same upstream channel to guarantee alignment. The analysis
    // output directory and the ROI YAML are staged as input files; the
    // notebook's parameters cell receives their staged names via params.yml.
    ch_report = ANALYSIS.out.outdir.combine(ch_roi_yaml.first()).map { meta, outdir, roi_yaml ->
        def parameters = [
            INDIR                  : outdir.name,
            SAMPLE_NAME            : meta.id,
            XENIUM_BUNDLE          : meta.samplesheet_xenium_bundle ?: '',
            SAMPLE_PUBLISHED_OUTDIR: meta.samplesheet_xenium_bundle ? "${params.outdir}/${params.mode}/qc/image_qc" : '',
            ROI_THRESHOLDS_YAML    : roi_yaml.name,
        ]
        [meta, parameters, [outdir, roi_yaml]]
    }
    REPORT(
        ch_report.map { meta, parameters, input_files -> [meta, notebook] },
        ch_report.map { meta, parameters, input_files -> parameters },
        ch_report.map { meta, parameters, input_files -> input_files },
        [],
    )

    emit:
    outdir = ANALYSIS.out.outdir // channel: [ val(meta), path(outdir) ]
    report = REPORT.out.html     // channel: [ val(meta), path(html) ]
}
