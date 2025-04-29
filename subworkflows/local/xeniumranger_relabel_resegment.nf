//
// run xeniumranger relabel & resegment to redine the xenium bundle
//

include { XENIUMRANGER_RELABEL   } from '../../modules/nf-core/xeniumranger/relabel/main'
include { XENIUMRANGER_RESEGMENT } from '../../modules/nf-core/xeniumranger/resegment/main'

workflow XENIUMRANGER_RELABEL_RESEGMENT {

    take:

    ch_bundle          // channel: [ val(meta), [ xenium-bundle-path ] ]
    ch_gene_panel      // channel: [ ["gene_panel.json"] ]

    main:

    ch_versions = Channel.empty()

    if ( params.relabel_genes ) {

        XENIUMRANGER_RELABEL ( ch_bundle, ch_gene_panel )
        ch_versions = ch_versions.mix ( XENIUMRANGER_RELABEL.out.versions )

        XENIUMRANGER_RESEGMENT ( XENIUMRANGER_RELABEL.out.bundle )
        ch_versions = ch_versions.mix ( XENIUMRANGER_RESEGMENT.out.versions )

    }


    emit:

    redefined_bundle = XENIUMRANGER_RESEGMENT.out.bundle // channel: [ val(meta), ["redefined-xenium-bundle"] ]

    versions = ch_versions                               // channel: [ versions.yml ]
}
