//
// Runs proseg for the xenium format and proseg2baysor to generate cell ploygons
//

include { PROSEG                           } from '../../../modules/local/proseg/preset/main'
include { PROSEG2BAYSOR                    } from '../../../modules/local/proseg/proseg2baysor/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/local/xeniumranger/import-segmentation/main'


workflow PROSEG_PRESET_PROSEG2BAYSOR {

    take:

    ch_bundle_path         // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    ch_transcripts_parquet // channel: [ val(meta), [ "transcripts.parquet" ] ]

    main:

    ch_versions = Channel.empty()

    // run proseg with the xenium format
    PROSEG ( ch_transcripts_parquet )
    ch_versions = ch_versions.mix( PROSEG.out.versions )

    // run proseg-to-baysor on the data generated with the proseg run
    PROSEG2BAYSOR ( 
        PROSEG.out.cell_polygons_2d.combine(PROSEG.out.transcript_metadata, by: 0)
    )
    ch_versions = ch_versions.mix( PROSEG2BAYSOR.out.versions )

    // run xeniumranger import-segmentation
    XENIUMRANGER_IMPORT_SEGMENTATION (
        ch_bundle_path
            .combine(PROSEG2BAYSOR.out.xr_polygons, by: 0)
            .combine(PROSEG2BAYSOR.out.xr_metadata, by: 0)
            .map {
                meta, bundle, xr_cell_polygons, xr_transcript_metadata -> tuple(
                    meta, // meta
                    bundle, // bundle
                    [], // coordinate_transform
                    [], // nuclei
                    [], // cells
                    xr_transcript_metadata, // transcript_assignment
                    xr_cell_polygons, // viz_polygons
                    "microns" // units
                )
            }

    )
    ch_versions = ch_versions.mix( XENIUMRANGER_IMPORT_SEGMENTATION.out.versions )

    emit:

    cell_polygons_2d      = PROSEG.out.cell_polygons_2d                 // channel: [ val(meta), [ "cell-polygons.geojson.gz" ] ]

    xr_polygons           = PROSEG2BAYSOR.out.xr_polygons               // channel: [ val(meta), [ "xr-cell-polygons.geojson" ] ]
    xr_metadata           = PROSEG2BAYSOR.out.xr_metadata               // channel: [ [ "xr-transcript-metadata.csv" ] ]
    coordinate_space      = ch_coordinate_space                         // channel: [ "microns" ]

    redefined_bundle      = XENIUMRANGER_IMPORT_SEGMENTATION.out.bundle // channel: [ val(meta), ["redefined-xenium-bundle"] ]

    versions              = ch_versions                                 // channel: [ versions.yml ]
}
