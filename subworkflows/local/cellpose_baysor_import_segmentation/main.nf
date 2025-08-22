//
// Run the cellpose, baysor and import-segmentation flow
//

include { RESOLIFT                         } from '../../../modules/local/resolift/main'
include { CELLPOSE as CELLPOSE_CELLS       } from '../../../modules/nf-core/cellpose/main'
include { CELLPOSE as CELLPOSE_NUCLEI      } from '../../../modules/nf-core/cellpose/main'
include { BAYSOR_RUN                       } from '../../../modules/local/baysor/run/main'
include { BAYSOR_PREPROCESS_TRANSCRIPTS    } from '../../../modules/local/baysor/preprocess/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow CELLPOSE_BAYSOR_IMPORT_SEGMENTATION {

    take:

    ch_morphology_image     // channel: [ val(meta), ["path-to-morphology.ome.tif"] ]
    ch_bundle_path          // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    ch_transcripts_parquet  // channel: [ val(meta), ["path-to-transcripts.parquet"] ]
    ch_config               // channel: ["path-to-xenium.toml"]

    main:

    ch_versions              = Channel.empty()
    ch_image                 = Channel.empty()
    ch_polygons              = Channel.empty()
    ch_segmentation          = Channel.empty()
    ch_filtered_transcripts  = Channel.empty()
    ch_cellpose_cells_mask   = Channel.empty()
    ch_cellpose_nuclei_mask  = Channel.empty()
    ch_cellpose_cells_cells  = Channel.empty()
    ch_cellpose_nuclei_cells = Channel.empty()
    ch_cellpose_cells_flows  = Channel.empty()
    ch_cellpose_nuclei_flows = Channel.empty()

    cellpose_model = params.cellpose_model ? (Channel.fromPath(params.cellpose_model, checkIfExists: true)) : []

    // sharpen morphology tiff if param - sharpen_tiff is true
    if ( params.sharpen_tiff ) {

        RESOLIFT ( ch_morphology_image )
        ch_versions = ch_versions.mix( RESOLIFT.out.versions )

        ch_image = RESOLIFT.out.enhanced_tiff

    } else {

        ch_image = ch_morphology_image

    }

    // run cellpose on the enhanced tiff
    if ( params.cell_segmentation_only ) {

        CELLPOSE_CELLS ( ch_image, cellpose_model, 'cells' )
        ch_versions = ch_versions.mix( CELLPOSE_CELLS.out.versions )

        ch_cellpose_cells_cells = CELLPOSE_CELLS.out.cells.map {
            _meta, cells -> return [ cells ]
        }
        ch_cellpose_cells_mask = CELLPOSE_CELLS.out.mask.map {
            _meta, mask -> return [ mask ]
        }
        ch_cellpose_cells_flows = CELLPOSE_CELLS.out.flows.map {
            _meta, flows -> return [ flows ]
        }

    }

    if ( params.nucleus_segmentation_only ) {

        CELLPOSE_NUCLEI ( ch_image, 'nuclei', 'nuclei' )
        ch_versions = ch_versions.mix( CELLPOSE_NUCLEI.out.versions )

        ch_cellpose_nuclei_cells = CELLPOSE_NUCLEI.out.cells.map {
            _meta, cells -> return [ cells ]
        }
        ch_cellpose_nuclei_mask = CELLPOSE_NUCLEI.out.mask.map {
            _meta, mask -> return [ mask ]
        }
        ch_cellpose_nuclei_flows = CELLPOSE_NUCLEI.out.flows.map {
            _meta, flows -> return [ flows ]
        }

    }


    // filter transcripts.parquet based on thresholds
    if ( params.filter_transcripts ) {

        BAYSOR_PREPROCESS_TRANSCRIPTS (
            ch_transcripts_parquet,
            params.min_qv,
            params.max_x,
            params.min_x,
            params.max_y,
            params.min_y
        )
        ch_versions = ch_versions.mix ( BAYSOR_PREPROCESS_TRANSCRIPTS.out.versions )

        ch_filtered_transcripts = BAYSOR_PREPROCESS_TRANSCRIPTS.out.transcripts_parquet

    }

    // run baysor with cellpose results
    if ( params.nucleus_segmentation_only ) {

        // run baysor with nuclei mask
        BAYSOR_RUN ( ch_filtered_transcripts, ch_cellpose_nuclei_mask, ch_config, 30 )
        ch_versions = ch_versions.mix ( BAYSOR_RUN.out.versions )

    } else if ( params.cell_segmentation_only ) {

        // run baysor with cell mask
        BAYSOR_RUN ( ch_filtered_transcripts, ch_cellpose_cells_mask, ch_config, 30 )
        ch_versions = ch_versions.mix ( BAYSOR_RUN.out.versions )

    } else {

        // run baysor with cell mask
        BAYSOR_RUN ( ch_filtered_transcripts, [], ch_config, 30 )
        ch_versions = ch_versions.mix ( BAYSOR_RUN.out.versions )

    }

    // run import-segmentation with baysor outs
    ch_segmentation = BAYSOR_RUN.out.segmentation.map {
        _meta, segmentation -> return [ segmentation ]
    }
    ch_polygons = BAYSOR_RUN.out.polygons2d


    XENIUMRANGER_IMPORT_SEGMENTATION (
        ch_bundle_path,
        [],
        [],
        [],
        ch_segmentation,
        ch_polygons,
        "microns"
    )
    ch_versions = ch_versions.mix ( XENIUMRANGER_IMPORT_SEGMENTATION.out.versions )

    emit:

    cells_mask  = ch_cellpose_cells_mask                           // channel: [ val(meta), [ "*masks.tif" ] ]
    cells_flows = ch_cellpose_cells_flows                          // channel: [ val(meta), [ "*flows.tif" ] ]
    cells_cells = ch_cellpose_cells_cells                          // channel: [ val(meta), [ "*seg.npy" ] ]
    nuclei_mask  = ch_cellpose_nuclei_mask                         // channel: [ val(meta), [ "*masks.tif" ] ]
    nuclei_flows = ch_cellpose_nuclei_flows                        // channel: [ val(meta), [ "*flows.tif" ] ]
    nuclei_cells = ch_cellpose_nuclei_cells                        // channel: [ val(meta), [ "*seg.npy" ] ]

    segmentation_mask = ch_segmentation                            // channel: [ val(meta), [ *segmentation.csv ] ]
    polygons          = ch_polygons                                // channel: [ val(meta), [ *polygons.json ] ]

    redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.bundle // channel: [ val(meta), ["redefined-xenium-bundle"] ]

    versions         = ch_versions                                 // channel: [ versions.yml ]
}

