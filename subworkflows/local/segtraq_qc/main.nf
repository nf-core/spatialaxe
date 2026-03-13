//
// Run SegTraQ QC
//

include { SEGTRAQ_BASELINE                     } from '../../../modules/local/segtraq/baseline/main'

workflow SEGTRAQ_QC {
    take:
    ch_spatialdata              // channel: [ val(meta), path("spatialdata.zarr") ]

    main:

    ch_versions = channel.empty()

    // run SegTraQ baseline QC metrics
    SEGTRAQ_BASELINE(ch_spatialdata)
    ch_versions = ch_versions.mix(SEGTRAQ_BASELINE.out.versions)

    emit:
    qc_results = SEGTRAQ_BASELINE.out.qc_results // channel: [ val(meta), path("segtraq_qc/*/") ]
    versions         = ch_versions // channel: [ versions.yml ]
}
