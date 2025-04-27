include { BIDCELL as BIDCELL_MODULE } from '../../modules/local/bidcell'

workflow BIDCELL {
    take:
    ch_samplesheet // This should contain: ID, bundle, image, ref, pos_markers, neg_markers

    main:
    ch_versions    = Channel.empty()

    ch_bidcell = ch_samplesheet.map {
        meta, bundle, image, ref, pos_markers, neg_markers ->
        [
            meta,
            bundle,
            image,
            // Optional. if missing, we provide a dummy file
            ref ?: file("${projectDir}/assets/bidcell_ref/sc_breast.csv"),
            pos_markers ?: file("${projectDir}/assets/bidcell_ref/sc_breast_markers_pos.csv"),
            neg_markers ?: file("${projectDir}/assets/bidcell_ref/sc_breast_markers_neg.csv"),
            ]
    }

    BIDCELL_MODULE(ch_bidcell)

    ch_versions = ch_versions.mix(BIDCELL_MODULE.out.versions)

    emit:
    outdir = BIDCELL_MODULE.out.outdir
    versions = ch_versions

}
