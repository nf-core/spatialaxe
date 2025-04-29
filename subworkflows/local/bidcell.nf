include { BIDCELL as BIDCELL_MODULE } from '../../modules/local/bidcell'

workflow BIDCELL {
    take:
    ch_samplesheet// This should contain: ID, bundle, image, ref, pos_markers, neg_markers

    main:
    ch_versions    = Channel.empty()

    ch_bidcell = ch_samplesheet.map {
        meta, bundle, image, ref, pos_markers, neg_markers ->
        [
            meta,
            // if parquet is not there, we assume csv.gz
            file(bundle + '/transcripts.parquet').exists() ? file(bundle + '/transcripts.parquet') : file(bundle + '/transcripts.csv.gz'),
            // if morphology focus is not there, we assume morphology mip
            file(bundle + '/morphology_focus/morphology_focus_0000.ome.tif').exists() ? file(bundle + '/morphology_focus/morphology_focus_0000.ome.tif') : file(bundle + '//morphology_mip.ome.tif'),
            // scRNA-seq ref set is optional. if missing, we provide a dummy one from BIDCELL github repo
            ref         ?: [],
            pos_markers ?: [],
            neg_markers ?: [],
        ]
    }

    BIDCELL_MODULE(ch_bidcell)

    ch_versions = ch_versions.mix(BIDCELL_MODULE.out.versions)

    emit:
    outdir = BIDCELL_MODULE.out.outdir
    versions = ch_versions

}
