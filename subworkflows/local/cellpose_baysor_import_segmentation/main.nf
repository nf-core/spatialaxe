//
// Run the cellpose, baysor and import-segmentation flow
//

include { RESOLIFT                         } from '../../../modules/local/resolift/main'
include { BAYSOR_RUN                       } from '../../../modules/local/baysor/run/main'
include { CELLPOSE as CELLPOSE_CELLS       } from '../../../modules/nf-core/cellpose/main'
include { EXTRACT_DAPI                     } from '../../../modules/local/utility/extract_dapi/main'
include { STARDIST as STARDIST_NUCLEI      } from '../../../modules/nf-core/stardist/main'
include { CONVERT_MASK_UINT32              } from '../../../modules/local/utility/convert_mask_uint32/main'
include { BAYSOR_PREPROCESS_TRANSCRIPTS    } from '../../../modules/local/baysor/preprocess/main'
include { RESIZE_TIF                       } from '../../../modules/local/utility/resize_tif/main'
include { GET_TRANSCRIPTS_COORDINATES      } from '../../../modules/local/utility/get_coordinates/main'
include { XENIUMRANGER_IMPORT_SEGMENTATION } from '../../../modules/nf-core/xeniumranger/import-segmentation/main'

workflow CELLPOSE_BAYSOR_IMPORT_SEGMENTATION {
    take:
    ch_morphology_image          // channel: [ val(meta), ["path-to-morphology.ome.tif"] ]
    ch_bundle_path               // channel: [ val(meta), ["path-to-xenium-bundle"] ]
    ch_transcripts_file       // channel: [ val(meta), ["path-to-transcripts.parquet"] ]
    ch_experiment_metadata       // channel: [ val(meta), ["path-to-experiment.xenium"] ]
    ch_config                    // channel: ["path-to-xenium.toml"]
    cell_segmentation_only       // value: bool
    cellpose_model               // value: path to cellpose model (or null)
    max_x                        // value: spatial filter upper x bound
    max_y                        // value: spatial filter upper y bound
    min_qv                       // value: minimum transcript QV
    min_x                        // value: spatial filter lower x bound
    min_y                        // value: spatial filter lower y bound
    nucleus_segmentation_only    // value: bool
    sharpen_tiff                 // value: bool
    stardist_nuclei_model        // value: stardist pretrained model name

    main:

    ch_transcripts = channel.empty()
    ch_imp_seg_inputs = channel.empty()
    ch_coordinate_space = channel.value("microns")


    // Use empty list when no model is provided; path input for official cellpose module
    cellpose_model_path = cellpose_model ? file(cellpose_model) : []
    stardist_model = stardist_nuclei_model ?: '2D_versatile_fluo'

    // sharpen morphology tiff if param - sharpen_tiff is true
    if (sharpen_tiff) {

        RESOLIFT(ch_morphology_image)

        ch_image = RESOLIFT.out.enhanced_tiff
    }
    else {

        ch_image = ch_morphology_image
    }


    // run cellpose on the morphology (enhanced) tiff
    if (cell_segmentation_only) {

        CELLPOSE_CELLS(ch_image, cellpose_model_path)
    }

    if (nucleus_segmentation_only) {

        // Extract DAPI channel, run StarDist, convert to uint32
        EXTRACT_DAPI(ch_image)

        STARDIST_NUCLEI(EXTRACT_DAPI.out.dapi, [stardist_model, []])

        CONVERT_MASK_UINT32(STARDIST_NUCLEI.out.mask)
    }


    // Always preprocess transcripts.parquet to CSV for Baysor 0.7.1 compatibility.
    // Baysor's Julia Parquet.jl cannot read zstd-compressed parquet files from Xenium bundles.
    // Also applies optional spatial/QV filtering when filter_transcripts is true.
    BAYSOR_PREPROCESS_TRANSCRIPTS(
        ch_transcripts_file,
        min_qv,
        max_x,
        min_x,
        max_y,
        min_y,
    )
    ch_transcripts = BAYSOR_PREPROCESS_TRANSCRIPTS.out.transcripts_file


    // run baysor with cellpose results
    if (nucleus_segmentation_only) {

        // check if the size of the segmentation mask matches the max transcripts coordinate range
        ch_resizetif_input = ch_transcripts
            .combine(CONVERT_MASK_UINT32.out.mask, by: 0)
            .combine(ch_experiment_metadata, by: 0)
            .map { meta, transcripts, mask, exp_meta ->
                tuple(
                    meta,
                    transcripts,
                    mask,
                    exp_meta,
                )
            }
        RESIZE_TIF(ch_resizetif_input)

        // run baysor with nuclei mask
        ch_baysor_input = ch_transcripts
            .combine(RESIZE_TIF.out.resized_mask, by: 0)
            .combine(ch_config)
            .map { meta, transcripts, mask, config ->
                tuple(
                    meta,
                    transcripts,
                    mask,
                    config,
                    30,
                )
            }
        BAYSOR_RUN(ch_baysor_input)
    }
    else if (cell_segmentation_only) {

        // check if the size of the segmentation mask matches the max transcripts coordinate range
        ch_resizetif_input = ch_transcripts
            .combine(CELLPOSE_CELLS.out.mask, by: 0)
            .combine(ch_experiment_metadata, by: 0)
            .map { meta, transcripts, mask, exp_meta ->
                tuple(
                    meta,
                    transcripts,
                    mask,
                    exp_meta,
                )
            }
        RESIZE_TIF(ch_resizetif_input)

        // run baysor with cell mask
        ch_baysor_input = ch_transcripts
            .combine(RESIZE_TIF.out.resized_mask, by: 0)
            .combine(ch_config)
            .map { meta, transcripts, mask, config ->
                tuple(
                    meta,
                    transcripts,
                    mask,
                    config,
                    30,
                )
            }
        BAYSOR_RUN(ch_baysor_input)
    }
    else {

        // run baysor without cell/nuclei mask
        ch_baysor_input = ch_transcripts
            .combine(ch_config)
            .map { meta, transcripts, config ->
                tuple(
                    meta,
                    transcripts,
                    [],
                    config,
                    30,
                )
            }
        BAYSOR_RUN(ch_baysor_input)
    }


    // run import-segmentation with baysor outs
    ch_imp_seg_inputs = ch_bundle_path
        .combine(BAYSOR_RUN.out.segmentation, by: 0)
        .map { meta, bundle, segmentation_csv, polygons2d ->
            tuple(
                meta,
                bundle,
                segmentation_csv,
                polygons2d,
                [],
                [],
                [],
                ch_coordinate_space.val,
            )
        }

    XENIUMRANGER_IMPORT_SEGMENTATION(ch_imp_seg_inputs)

    emit:
    coordinate_space = ch_coordinate_space                         // channel: [ val("microns") ]
    redefined_bundle = XENIUMRANGER_IMPORT_SEGMENTATION.out.outs // channel: [ val(meta), ["redefined-xenium-bundle"] ]
}
