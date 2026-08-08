//
// TRANSCRIPT_QC: transcript / molecule-level quality control for a Xenium bundle.
//
// Runs the transcript QC analysis (per-transcript and per-cell QC metrics +
// figures) and renders an HTML report from the analysis outputs with the
// nf-core QUARTONOTEBOOK module. Versions are reported via the `versions`
// topic channel by each module.
//

include { TRANSCRIPT_QC_PROCESSING as ANALYSIS } from '../../../modules/local/transcript_qc/main'
include { QUARTONOTEBOOK as REPORT             } from '../../../modules/nf-core/quartonotebook/main'

workflow TRANSCRIPT_QC {
    take:
    ch_input // channel: [ val(meta), val(parameters), path(input_files) ]
    notebook // path: the transcript QC report .qmd

    main:
    // Run computational processing
    ANALYSIS(ch_input)

    // Render the report with QUARTONOTEBOOK. Its four inputs are separate
    // channels paired by emission order, so all per-sample channels are derived
    // from the same upstream channel to guarantee alignment. The analysis
    // output directory is staged as an input file; the notebook's parameters
    // cell receives its staged name via params.yml.
    ch_report = ANALYSIS.out.outdir.map { meta, outdir ->
        def parameters = [
            INDIR                  : outdir.name,
            SAMPLE_NAME            : meta.id,
            XENIUM_BUNDLE          : meta.samplesheet_xenium_bundle ?: '',
            SAMPLE_PUBLISHED_OUTDIR: meta.samplesheet_xenium_bundle ? "${params.outdir}/${params.mode}/qc/transcript_qc" : '',
        ]
        [meta, parameters, outdir]
    }
    REPORT(
        ch_report.map { meta, parameters, outdir -> [meta, notebook] },
        ch_report.map { meta, parameters, outdir -> parameters },
        ch_report.map { meta, parameters, outdir -> outdir },
        [],
    )

    emit:
    outdir = ANALYSIS.out.outdir // channel: [ val(meta), path(outdir) ]
    report = REPORT.out.html     // channel: [ val(meta), path(html) ]
}
