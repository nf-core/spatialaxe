//
// TRANSCRIPT_QC: transcript / molecule-level quality control for a Xenium bundle.
//
// Runs the transcript QC analysis (per-transcript and per-cell QC metrics +
// figures) and renders an HTML report from the analysis outputs with Quarto.
// Versions are reported via the `versions` topic channel by each module.
//

include { TRANSCRIPT_QC_PROCESSING as ANALYSIS } from '../../../modules/local/transcript_qc/main'
include { QUARTO as REPORT } from '../../../modules/local/quarto/main'

workflow TRANSCRIPT_QC {
    take:
    ch_input    // channel: [ val(meta), val(parameters), path(input_files) ]
    ch_notebook // channel: [ val(meta2), path(qmd_file) ]

    main:
    // Run computational processing
    ANALYSIS(ch_input)

    // Generate report using the Quarto module. Transcript QC has no ROI
    // thresholds, so the ROI YAML input is an empty placeholder.
    REPORT(
        ANALYSIS.out.outdir,
        ch_notebook,
        [],
    )

    emit:
    outdir = ANALYSIS.out.outdir // channel: [ val(meta), path(outdir) ]
    report = REPORT.out.report   // channel: [ val(meta), path(html) ]
}
