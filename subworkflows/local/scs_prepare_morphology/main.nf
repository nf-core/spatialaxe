//
// Prepare SCS-compatible input from Xenium transcripts and pass morphology for downstream SCS segmentation
//

include { XENIUM2SCS } from '../../../modules/local/utility/xenium2scs/main'

workflow SCS_PREPARE_MORPHOLOGY {
    take:
    ch_morphology_image    // channel: [ val(meta), ["path-to-morphology.ome.tif"] ]
    ch_transcripts_parquet // channel: [ val(meta), ["path-to-transcripts.parquet"] ]
    ch_experiment_xenium   // channel: [ val(meta), ["path-to-experiment.xenium"] ]

    main:

    ch_versions = channel.empty()

    // convert Xenium transcripts.parquet to SCS tabular input format
    ch_xenium2scs_input = ch_transcripts_parquet
        .join(ch_morphology_image, by: 0)
        .join(ch_experiment_xenium, by: 0)
    XENIUM2SCS(ch_xenium2scs_input)
    ch_versions = ch_versions.mix(XENIUM2SCS.out.versions)

    emit:
    morphology_image = ch_morphology_image      // channel: [ val(meta), ["path-to-morphology.ome.tif"] ]
    morphology_2d    = XENIUM2SCS.out.morph2d_tif // channel: [ val(meta), ["path-to-morph2d.tif"] ]
    scs_input_bgi_tsv = XENIUM2SCS.out.scs_input_bgi_tsv // channel: [ val(meta), ["path-to-scs_input_bgi.tsv"] ]
    metrics          = XENIUM2SCS.out.metrics   // channel: [ val(meta), ["path-to-xenium2scs_metrics.tsv"] ]
    versions         = ch_versions              // channel: [ versions.yml ]
}
