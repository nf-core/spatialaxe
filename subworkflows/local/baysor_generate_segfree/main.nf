//
// Run baysor segfree
//

include { BAYSOR_SEGFREE } from '../../../modules/local/baysor/segfree/main'


workflow BAYSOR_GENERATE_SEGFREE {

    take:

    ch_transcripts_parquet // channel: [ val(meta), ["transcripts.parquet"] ]
    ch_config              // channel: [ ["path-to-xenium.toml"] ]

    main:

    ch_versions = Channel.empty()

    // run baysor segfree
    BAYSOR_SEGFREE (
        ch_transcripts_parquet,
        ch_config
    )
    ch_versions = ch_versions.mix ( BAYSOR_SEGFREE.out.versions )

    emit:

    ncvs     = BAYSOR_SEGFREE.out.ncvs // channel: [ val(meta), ["ncvs.loom"] ]

    versions = ch_versions             // channel: [ versions.yml ]
}
