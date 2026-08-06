//
// IMAGE_QC: image-based quality control for a Xenium bundle.
//
// Runs the image QC analysis (focus / SNR / morphology metrics + figures) and
// renders an HTML report from the analysis outputs with Quarto. Versions are
// reported via the `versions` topic channel by each module (collected centrally
// by the pipeline), so this subworkflow does not thread versions through emit.
//

include { IMAGE_QC_ANALYSIS as ANALYSIS } from '../../../modules/local/image_qc/main'
include { QUARTO as REPORT } from '../../../modules/local/quarto/main'

workflow IMAGE_QC {
    take:
    ch_input    // channel: [ val(meta), val(parameters), path(input_files) ]
    ch_notebook // channel: [ val(meta2), path(qmd_file) ]

    main:
    // ROI threshold config: user-provided path, else the bundled default.
    ch_roi_yaml = channel.fromPath(
        params.roi_image_qc_thresholds_yaml ?: "${projectDir}/conf/roi_image_qc_thresholds.yaml",
        checkIfExists: true,
    )

    // Single process handles both ROI-based and cell-level QC
    ANALYSIS(ch_input, ch_roi_yaml.first())

    // Generate report using the Quarto module
    REPORT(
        ANALYSIS.out.outdir,
        ch_notebook,
        ch_roi_yaml.first(),
    )

    emit:
    outdir = ANALYSIS.out.outdir // channel: [ val(meta), path(outdir) ]
    report = REPORT.out.report   // channel: [ val(meta), path(html) ]
}
