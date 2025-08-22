//
// Runs proseg for the xenium format and proseg2baysor to generate cell ploygons
//

include { PROSEG                           } from '../../../modules/local/proseg/preset/main'
include { PROSEG2BAYSOR                    } from '../../../modules/local/proseg/proseg2baysor/main'
include { PARQUET_TO_CSV                   } from '../../../modules/local/spatialconverter/parquet_to_csv/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow PROSEG_PRESET_PROSEG2BAYSOR {

    take:

    ch_bundle_path         // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    ch_transcripts_parquet // channel: [ val(meta), [ "transcripts.parquet" ] ]

    main:

    ch_versions = Channel.empty()

    // run parquet-to-csv
    PARQUET_TO_CSV ( ch_transcripts_parquet, ".gz" )
    ch_versions = ch_versions.mix( PARQUET_TO_CSV.out.versions )

    // run proseg with the xenium format
    PROSEG ( PARQUET_TO_CSV.out.transcripts_csv )
    ch_versions = ch_versions.mix( PROSEG.out.versions )

    // run proseg-to-baysor on the data generated with the proseg run
    PROSEG2BAYSOR ( PROSEG.out.cell_polygons_2d, PROSEG.out.transcript_metadata )
    ch_versions = ch_versions.mix( PROSEG2BAYSOR.out.versions )

    ch_metadata = PROSEG2BAYSOR.out.xr_metadata
    ch_polygons = PROSEG2BAYSOR.out.xr_polygons.map {
        _meta, polygons -> return [ polygons ]
    }

    // run xeniumranger import-segmentation
    XENIUMRANGER_IMPORT_SEGMENTATION (
        ch_bundle_path,
        [],
        [],
        [],
        ch_metadata,
        ch_polygons,
        "microns"
    )
    ch_versions = ch_versions.mix( XENIUMRANGER_IMPORT_SEGMENTATION.out.versions )

    emit:

    cell_polygons_2d      = PROSEG.out.cell_polygons_2d                 // channel: [ val(meta), [ "cell-polygons.geojson.gz" ] ]

    xr_polygons           = PROSEG2BAYSOR.out.xr_polygons               // channel: [ val(meta), [ "xr-cell-polygons.geojson" ] ]
    xr_metadata           = PROSEG2BAYSOR.out.xr_metadata               // channel: [ [ "xr-transcript-metadata.csv" ] ]

    redefined_bundle      = XENIUMRANGER_IMPORT_SEGMENTATION.out.bundle // channel: [ val(meta), ["redefined-xenium-bundle"] ]

    versions              = ch_versions                                 // channel: [ versions.yml ]
}

