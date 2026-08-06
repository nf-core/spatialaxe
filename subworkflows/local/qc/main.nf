//
// QC: quality-control layer for a Xenium bundle.
//
// Thin wrapper that runs the two QC subworkflows — image QC and transcript QC —
// on a reported bundle and emits their analysis directories and HTML reports.
// It deliberately holds no pre/post-segmentation logic and no MultiQC: when and
// on which bundle QC runs is decided by the caller (workflows/spatialaxe.nf).
//

include { IMAGE_QC } from '../image_qc/main'
include { TRANSCRIPT_QC } from '../transcript_qc/main'

workflow QC {
    take:
    ch_bundle // channel: [ val(meta), path(xenium_bundle) ]

    main:
    // Shared per-sample meta for both QC reports.
    ch_qc_meta = ch_bundle.map { meta, bundle ->
        [
            [
                id: meta.id,
                samplesheet_xenium_bundle: bundle.toString(),
                xenium_bundle_source: meta.xenium_bundle_source ?: '',
                cropped: meta.cropped ?: false,
            ],
            bundle,
        ]
    }

    // ---------------------------- Image QC ------------------------------
    imageqc_input_ch = ch_qc_meta.map { meta, bundle ->
        tuple(
            meta,
            // parameters (flat key/value list, collated by the module)
            ["STAIN_NAMES", params.stain_names],
            // input files
            [bundle],
        )
    }
    image_qc_notebook = file("${projectDir}/assets/notebooks/xenium_image_qc_report.qmd", checkIfExists: true)
    IMAGE_QC(
        imageqc_input_ch,
        [[id: image_qc_notebook.getBaseName()], image_qc_notebook],
    )

    // -------------------------- Transcript QC ---------------------------
    transcriptqc_input_ch = ch_qc_meta.map { meta, bundle ->
        tuple(
            meta,
            // parameters (NON_GENE_PREFIX is the key the module's arg builder reads)
            ["NON_GENE_PREFIX", params.neg_control_prefix],
            // input files
            [bundle],
        )
    }
    transcript_qc_notebook = file("${projectDir}/assets/notebooks/transcript_qc.qmd", checkIfExists: true)
    TRANSCRIPT_QC(
        transcriptqc_input_ch,
        [[id: transcript_qc_notebook.getBaseName()], transcript_qc_notebook],
    )

    emit:
    image_qc_outdir      = IMAGE_QC.out.outdir       // channel: [ val(meta), path(outdir) ]
    image_qc_report      = IMAGE_QC.out.report       // channel: [ val(meta), path(html) ]
    transcript_qc_outdir = TRANSCRIPT_QC.out.outdir  // channel: [ val(meta), path(outdir) ]
    transcript_qc_report = TRANSCRIPT_QC.out.report  // channel: [ val(meta), path(html) ]
}
