//
// generate spatialdata object from the spatialxe layers
//

include { SPATIALDATA_WRITE as SPATIALDATA_WRITE_RAW_BUNDLE       } from '../../../modules/local/spatialdata/write/main'
// spoQC general stuff
include { SPOQC_ANNOTATION       } from '../../../modules/local/spoQC/annotation/main'
include { SPOQC_WHOLE_SLIDE       } from '../../../modules/local/spoQC/whole_slide/main'
include { SPOQC_GENERAL       } from '../../../modules/local/spoQC/general/main'
include { SPOQC_BUBBLE       } from '../../../modules/local/spoQC/bubble/main'
include { SPOQC_DOUBLET       } from '../../../modules/local/spoQC/doublet/main'
include { SPOQC_VOID       } from '../../../modules/local/spoQC/void/main'
include { SPOQC_CELL       } from '../../../modules/local/spoQC/cell/main'
include { SPOQC_AMBIENT    } from '../../../modules/local/spoQC/ambient/main'
// spoQC hqcr
include { SPOQC_HQCR_IDENT       } from '../../../modules/local/spoQC/hqcr_ident/main'
include { SPOQC_HQCR_CELLTYPE       } from '../../../modules/local/spoQC/hqcr_celltype/main'
// spoQC hqpr
include { SPOQC_HQPR_METRICES       } from '../../../modules/local/spoQC/hqpr_metrices/main'
include { SPOQC_HQPR_CLUSTERING       } from '../../../modules/local/spoQC/hqpr_clustering/main'
include { SPOQC_HQPR_REFINEMENT       } from '../../../modules/local/spoQC/hqpr_refinement/main'
include { SPOQC_HQPR_BOUNDING_BOX       } from '../../../modules/local/spoQC/hqpr_bounding_box/main'
include { SPOQC_HQPR_CELLTYPE       } from '../../../modules/local/spoQC/hqpr_celltype/main'
// spoQC hqtr
include { SPOQC_HQTR_METRICES       } from '../../../modules/local/spoQC/hqtr_metrices/main'
include { SPOQC_HQTR_AC       } from '../../../modules/local/spoQC/hqtr_ac/main'
include { SPOQC_HQTR_QV       } from '../../../modules/local/spoQC/hqtr_qv/main'
include { SPOQC_HQTR_CLUSTERING       } from '../../../modules/local/spoQC/hqtr_clustering/main'
include { SPOQC_HQTR_REFINEMENT       } from '../../../modules/local/spoQC/hqtr_refinement/main'
include { SPOQC_HQTR_BOUNDING_BOX       } from '../../../modules/local/spoQC/hqtr_bounding_box/main'
include { SPOQC_HQTR_CELLTYPE       } from '../../../modules/local/spoQC/hqtr_celltype/main'
// spoQC downstream
include { SPOQC_COMBINE_MASKS       } from '../../../modules/local/spoQC/combine_masks/main'
include { SPOQC_TRANSCRIPT       } from '../../../modules/local/spoQC/transcript/main'
include { SPOQC_CELLCYCLE       } from '../../../modules/local/spoQC/cellcycle/main'
include { SPOQC_MODEL       } from '../../../modules/local/spoQC/model/main'
include { SPOQC_MARKER       } from '../../../modules/local/spoQC/marker/main'
// spoQC final analysis
include { SPOQC_ANALYSIS_OVERVIEW       } from '../../../modules/local/spoQC/analysis_overview/main'
include { SPOQC_ANALYSIS_CATEGORY       } from '../../../modules/local/spoQC/analysis_category/main'
include { SPOQC_ANALYSIS_CLUSTER       } from '../../../modules/local/spoQC/analysis_cluster/main'
// spoQC final report
include { SPOQC_FINALREPORT       } from '../../../modules/local/spoQC/finalreport/main'

workflow SPOQC {

    take:
    ch_bundle_path          // channel: [ val(meta), [ "path-to-xenium-bundle" ] ]
    ch_annotation_src       // channel: [ [ "path-to-annotation-file" ] ]
    ch_stainings            // channel: [ [ 1,2,.... ] ]

    main:

    if ( params.nucleus_segmentation_only && params.cell_segmentation_only ) {

        ch_segmented_object = channel.value('cells_and_nuclei')

    }       

    else if ( params.nucleus_segmentation_only ) {

        ch_segmented_object = channel.value('nuclei')

    }

    else if ( params.cell_segmentation_only ) {

        ch_segmented_object = channel.value('cells')

    } else {

        ch_segmented_object = channel.value([])

    }

    ch_coordinate_space = channel.value("all")

    // write spatialdata object from the raw xenium bundle
    SPATIALDATA_WRITE_RAW_BUNDLE (
        ch_bundle_path,
        'spatialdata_raw',
        ch_segmented_object,
        ch_coordinate_space,
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // General
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // keep only actual, usable annotations (non-null and not an empty list)
    ch_annotation_present = ch_annotation_src.filter { a -> a && (!(a instanceof List) || !a.isEmpty()) }

    SPOQC_ANNOTATION(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        ch_annotation_src,
        "annotation",
    )
    ch_annotation_path = ch_annotation_present.mix( SPOQC_ANNOTATION.out.annotation )

    SPOQC_GENERAL(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        ch_annotation_path,
        "generalqc",
    )

    SPOQC_WHOLE_SLIDE(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "whole_slide_qc",
    )

    SPOQC_BUBBLE(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "bubbleqc",
    )

    SPOQC_DOUBLET(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        ch_annotation_path,
        "doubletqc",
    )

    SPOQC_VOID(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "voidqc",
    )

    SPOQC_CELL(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "cellqc",
    )

    SPOQC_AMBIENT(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "ambientqc",
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HQCR
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    SPOQC_HQCR_IDENT(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "hqcr_ident",
        SPOQC_GENERAL.out.tmp,
        SPOQC_BUBBLE.out.tmp,
        SPOQC_DOUBLET.out.tmp,
        SPOQC_VOID.out.tmp,
        SPOQC_CELL.out.tmp,
    )

    SPOQC_HQCR_CELLTYPE(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        ch_annotation_path,
        "hqcr_celltype",
        SPOQC_GENERAL.out.tmp,
        SPOQC_BUBBLE.out.tmp,
        SPOQC_DOUBLET.out.tmp,
        SPOQC_VOID.out.tmp,
        SPOQC_CELL.out.tmp,
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HQPR
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    SPOQC_HQPR_METRICES(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        ch_stainings,
        "hqpr_metrices",
    )

    ch_spatialdata_stainings = SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata.combine(ch_stainings)
    ch_annotation_stainings = ch_annotation_path.combine(ch_stainings)

    SPOQC_HQPR_CLUSTERING(
        ch_spatialdata_stainings,
        "hqpr_clustering",
        SPOQC_HQPR_METRICES.out.metrices,
    )

    SPOQC_HQPR_REFINEMENT(
        ch_spatialdata_stainings,
        "hqpr_refinement",  
        SPOQC_HQPR_CLUSTERING.out.mask,
    )

    SPOQC_HQPR_BOUNDING_BOX(
        ch_spatialdata_stainings,
        "hqpr_bounding_box",  
        SPOQC_HQPR_REFINEMENT.out.mask_smoothed,
    )

    // ch_masks_joined emits: tuple(val(staining), path(mask), path(mask))
    // ch_masks_joined = SPOQC_HQPR_CLUSTERING.out.mask.join( SPOQC_HQPR_REFINEMENT.out.mask_smoothed )

    // SPOQC_HQPR_CELLTYPE(
    //     ch_spatialdata_stainings,
    //     ch_annotation_stainings,
    //     "hqpr_celltype",  
    //     SPOQC_GENERAL.out.tmp.combine(ch_stainings),
    //     SPOQC_BUBBLE.out.tmp.combine(ch_stainings),
    //     SPOQC_DOUBLET.out.tmp.combine(ch_stainings),
    //     SPOQC_VOID.out.tmp.combine(ch_stainings),
    //     SPOQC_CELL.out.tmp.combine(ch_stainings),
    //     ch_masks_joined,
    // )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // HQTR
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPOQC_HQTR_METRICES(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "hqtr_metrices",
    )

    SPOQC_HQTR_AC(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "hqtr_ac",
        SPOQC_AMBIENT.out.tmp,
    )

    SPOQC_HQTR_QV(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "hqtr_qv",
    )

    SPOQC_HQTR_CLUSTERING(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "hqtr_clustering",
        SPOQC_HQTR_METRICES.out.metrices,
        SPOQC_HQTR_QV.out.tmp,
        SPOQC_HQTR_AC.out.tmp,
    )

    SPOQC_HQTR_REFINEMENT(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "hqtr_refinement",  
        SPOQC_HQTR_CLUSTERING.out.mask,
    )

    SPOQC_HQTR_BOUNDING_BOX(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "hqtr_bounding_box",  
        SPOQC_HQTR_REFINEMENT.out.mask_smoothed,
    )

    // SPOQC_HQTR_CELLTYPE(
    //     SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
    //     ch_annotation_path,
    //     "hqtr_celltype",  
    //     SPOQC_GENERAL.out.tmp,
    //     SPOQC_BUBBLE.out.tmp,
    //     SPOQC_DOUBLET.out.tmp,
    //     SPOQC_VOID.out.tmp,
    //     SPOQC_CELL.out.tmp,
    //     SPOQC_HQTR_CLUSTERING.out.mask,
    //     SPOQC_HQTR_REFINEMENT.out.mask_smoothed,
    // )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Downstream
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    SPOQC_COMBINE_MASKS(
        ch_spatialdata_stainings,
        "combine_masks",
        SPOQC_HQCR_IDENT.out.mask.combine(ch_stainings),
        SPOQC_HQPR_CLUSTERING.out.mask,
        SPOQC_HQTR_CLUSTERING.out.mask.combine(ch_stainings),
        SPOQC_HQCR_IDENT.out.mask_smoothed.combine(ch_stainings),
        SPOQC_HQPR_REFINEMENT.out.mask_smoothed,
        SPOQC_HQTR_REFINEMENT.out.mask_smoothed.combine(ch_stainings),
    )

    SPOQC_TRANSCRIPT(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        ch_annotation_path,
        "transcriptqc",
    )

    SPOQC_CELLCYCLE(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "cellcycleqc",
    )

    SPOQC_MODEL(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "modelqc",
    )

    // SPOQC_MARKER(
    //     SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
    //     ch_annotation_path,
    //     "markerqc"
    // )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Analysis
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ch_files_hqpr_metrics = SPOQC_HQPR_METRICES.out.metrices
        .map { _idx, p -> p }
        .collect()

    ch_files_hqpr_masks = SPOQC_HQPR_CLUSTERING.out.mask
        .map { _idx, p -> p }
        .collect()

    ch_files_hqpr_masks_smoothed = SPOQC_HQPR_REFINEMENT.out.mask_smoothed
        .map { _idx, p -> p }
        .collect()

    SPOQC_ANALYSIS_OVERVIEW(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        ch_annotation_path,
        "analysis_overview",
        SPOQC_GENERAL.out.tmp,
        SPOQC_BUBBLE.out.tmp,
        SPOQC_DOUBLET.out.tmp,
        SPOQC_VOID.out.tmp,
        SPOQC_CELL.out.tmp,
        SPOQC_HQCR_IDENT.out.mask,
        SPOQC_HQCR_IDENT.out.mask_smoothed,
        SPOQC_HQTR_QV.out.tmp,
        SPOQC_HQTR_AC.out.tmp,
        SPOQC_HQTR_METRICES.out.metrices,
        SPOQC_HQTR_REFINEMENT.out.mask_smoothed,
        SPOQC_HQTR_CLUSTERING.out.mask,
        ch_files_hqpr_metrics,
        ch_files_hqpr_masks_smoothed,
        ch_files_hqpr_masks,
    )

    SPOQC_ANALYSIS_CATEGORY(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        ch_annotation_path,
        "analysis_category",
        SPOQC_GENERAL.out.tmp,
        SPOQC_BUBBLE.out.tmp,
        SPOQC_DOUBLET.out.tmp,
        SPOQC_VOID.out.tmp,
        SPOQC_CELL.out.tmp,
        SPOQC_HQCR_IDENT.out.mask,
        SPOQC_HQCR_IDENT.out.mask_smoothed,
        SPOQC_HQTR_QV.out.tmp,
        SPOQC_HQTR_AC.out.tmp,
        SPOQC_HQTR_METRICES.out.metrices,
        SPOQC_HQTR_REFINEMENT.out.mask_smoothed,
        SPOQC_HQTR_CLUSTERING.out.mask,
        ch_files_hqpr_metrics,
        ch_files_hqpr_masks_smoothed,
        ch_files_hqpr_masks,
    )

    SPOQC_ANALYSIS_CLUSTER(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        ch_annotation_path,
        "analysis_cluster",
        SPOQC_GENERAL.out.tmp,
        SPOQC_BUBBLE.out.tmp,
        SPOQC_DOUBLET.out.tmp,
        SPOQC_VOID.out.tmp,
        SPOQC_CELL.out.tmp,
        SPOQC_HQCR_IDENT.out.mask,
        SPOQC_HQCR_IDENT.out.mask_smoothed,
        SPOQC_HQTR_QV.out.tmp,
        SPOQC_HQTR_AC.out.tmp,
        SPOQC_HQTR_METRICES.out.metrices,
        SPOQC_HQTR_REFINEMENT.out.mask_smoothed,
        SPOQC_HQTR_CLUSTERING.out.mask,
        ch_files_hqpr_metrics,
        ch_files_hqpr_masks_smoothed,
        ch_files_hqpr_masks,
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Final Report
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // collect the per-staining report fragments into a single list per sample
    ch_report_hqpr_metrices     = SPOQC_HQPR_METRICES.out.report.map { _staining, p -> p }.collect()
    ch_report_hqpr_clustering   = SPOQC_HQPR_CLUSTERING.out.report.map { _staining, p -> p }.collect()
    ch_report_hqpr_refinement   = SPOQC_HQPR_REFINEMENT.out.report.map { _staining, p -> p }.collect()
    ch_report_hqpr_bounding_box = SPOQC_HQPR_BOUNDING_BOX.out.report.map { _staining, p -> p }.collect()
    ch_report_combine_masks     = SPOQC_COMBINE_MASKS.out.report.collect()

    SPOQC_FINALREPORT(
        SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata,
        "final_report",
        SPOQC_GENERAL.out.report,
        SPOQC_DOUBLET.out.report,
        SPOQC_VOID.out.report,
        SPOQC_CELL.out.report,
        SPOQC_HQCR_IDENT.out.report,
        SPOQC_HQCR_CELLTYPE.out.report,
        ch_report_hqpr_metrices,
        ch_report_hqpr_clustering,
        ch_report_hqpr_refinement,
        ch_report_hqpr_bounding_box,
        // SPOQC_HQPR_CELLTYPE.out.report,
        SPOQC_HQTR_METRICES.out.report,
        SPOQC_HQTR_AC.out.report,
        SPOQC_HQTR_QV.out.report,
        SPOQC_HQTR_CLUSTERING.out.report,
        SPOQC_HQTR_REFINEMENT.out.report,
        SPOQC_HQTR_BOUNDING_BOX.out.report,
        // SPOQC_HQTR_CELLTYPE.out.report,
        ch_report_combine_masks,
        SPOQC_TRANSCRIPT.out.report,
        SPOQC_CELLCYCLE.out.report,
        SPOQC_MODEL.out.report,
        SPOQC_ANALYSIS_OVERVIEW.out.report,
        SPOQC_ANALYSIS_CATEGORY.out.report,
        SPOQC_ANALYSIS_CLUSTER.out.report,
    )

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Output
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    emit:

    ch_sd_raw       = SPATIALDATA_WRITE_RAW_BUNDLE.out.spatialdata         // channel: [ val(meta), "spatialdata_raw" ]
}
