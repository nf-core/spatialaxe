//
// run xeniumranger relabel & resegment to redine the xenium bundle
//

include { XENIUMRANGER_RELABEL   } from '../../../modules/nf-core/xeniumranger/relabel/main'
include { XENIUMRANGER_RESEGMENT } from '../../../modules/nf-core/xeniumranger/resegment/main'

workflow XENIUMRANGER_RELABEL_RESEGMENT {
    take:
    ch_bundle_path // channel: [ val(meta), [ "path-to-xenium-bundle" ] ]
    ch_gene_panel  // channel: [ val(meta), ["path-to-gene_panel.json"] ]

    main:

    ch_versions = Channel.empty()

    // Combine bundle path with gene panel into a single tuple for relabel
    XENIUMRANGER_RELABEL(
        ch_bundle_path.combine(ch_gene_panel, by: 0),
    )

    XENIUMRANGER_RESEGMENT(
        XENIUMRANGER_RELABEL.out.outs
    )

    emit:
    redefined_bundle = XENIUMRANGER_RESEGMENT.out.outs // channel: [ val(meta), ["redefined-xenium-bundle"] ]
    versions         = ch_versions // channel: [ versions.yml ]
}
