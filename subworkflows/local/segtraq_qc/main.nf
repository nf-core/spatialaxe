//
// Run SegTraQ QC
//

include { SEGTRAQ_BASELINE                     } from '../../../modules/local/segtraq/baseline/main'
include { SEGTRAQ_CLUSTERING_STABILITY         } from '../../../modules/local/segtraq/clustering_stability/main'
include { SEGTRAQ_REGION_SIMILARITY            } from '../../../modules/local/segtraq/region_similarity/main'
include { SEGTRAQ_VOLUME                       } from '../../../modules/local/segtraq/volume/main'

workflow SEGTRAQ_QC {
    take:
    ch_spatialdata              // channel: [ val(meta), path("spatialdata.zarr") ]

    main:

    ch_versions = channel.empty()
    ch_qc = channel.empty()
    def modules_to_run = params.segtraq_modules == 'all' ?
    ['baseline', 'clustering_stability', 'region_similarity', 'volume'] : params.segtraq_modules.tokenize(',')


    // run SegTraQ baseline QC metrics
    if ('baseline' in modules_to_run) {
        SEGTRAQ_BASELINE(ch_spatialdata)
        ch_versions = ch_versions.mix(SEGTRAQ_BASELINE.out.versions)
        ch_qc = ch_qc.mix(SEGTRAQ_BASELINE.out.qc_results)}

    if ('clustering_stability' in modules_to_run) {
        SEGTRAQ_CLUSTERING_STABILITY(ch_spatialdata)
        ch_versions = ch_versions.mix(SEGTRAQ_CLUSTERING_STABILITY.out.versions)
        ch_qc = ch_qc.mix(SEGTRAQ_CLUSTERING_STABILITY.out.qc_results)}

    if ('region_similarity' in modules_to_run) {
        SEGTRAQ_REGION_SIMILARITY(ch_spatialdata)
        ch_versions = ch_versions.mix(SEGTRAQ_REGION_SIMILARITY.out.versions)
        ch_qc = ch_qc.mix(SEGTRAQ_REGION_SIMILARITY.out.qc_results)}

    if ('volume' in modules_to_run) {
        SEGTRAQ_VOLUME(ch_spatialdata)
        ch_versions = ch_versions.mix(SEGTRAQ_VOLUME.out.versions)
        ch_qc = ch_qc.mix(SEGTRAQ_VOLUME.out.qc_results)
    }

    emit:
    qc_results = ch_qc // channel: [ val(meta), path("segtraq_qc/*/") ]
    versions         = ch_versions // channel: [ versions.yml ]
}
