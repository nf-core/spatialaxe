//
// Run SegTraQ QC
//

include { SEGTRAQ_BASELINE                     } from '../../../modules/local/segtraq/baseline/main'
include { SEGTRAQ_CLUSTERING_STABILITY         } from '../../../modules/local/segtraq/clustering_stability/main'

workflow SEGTRAQ_QC {
    take:
    ch_spatialdata              // channel: [ val(meta), path("spatialdata.zarr") ]

    main:

    ch_versions = channel.empty()
    ch_qc = channel.empty()

    // run SegTraQ baseline QC metrics
    SEGTRAQ_BASELINE(ch_spatialdata)
    ch_versions = ch_versions.mix(SEGTRAQ_BASELINE.out.versions)
    ch_qc = ch_qc.mix(SEGTRAQ_BASELINE.out.qc_results)

    SEGTRAQ_CLUSTERING_STABILITY(ch_spatialdata)
    ch_versions = ch_versions.mix(SEGTRAQ_CLUSTERING_STABILITY.out.versions)
    ch_qc = ch_qc.mix(SEGTRAQ_CLUSTERING_STABILITY.out.qc_results)

    emit:
    qc_results = ch_qc // channel: [ val(meta), path("segtraq_qc/*/") ]
    versions         = ch_versions // channel: [ versions.yml ]
}
