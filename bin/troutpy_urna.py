#!/usr/bin/env python3
"""Run troutpy uRNA analysis on SpatialData (Basic tutorial workflow).

Called from the TROUTPY Nextflow module. Must not be named troutpy.py — that
shadows the troutpy package on sys.path when executed from bin/.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

DEFAULT_CONTROL_CODEWORDS = [
    "negative_control_probe",
    "unassigned_codeword",
    "deprecated_codeword",
    "genomic_control_probe",
    "negative_control_codeword",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Input SpatialData Zarr")
    parser.add_argument("--output-zarr", type=Path, required=True, help="Output SpatialData Zarr")
    parser.add_argument(
        "--plots-dir",
        type=Path,
        default=Path("report/troutpy/plots"),
        help="Directory for saved figures",
    )
    parser.add_argument(
        "--tables-dir",
        type=Path,
        default=Path("report/troutpy/tables"),
        help="Directory for exported tables (CSV)",
    )
    parser.add_argument("--skip-plots", action="store_true", help="Skip plotting")
    parser.add_argument(
        "--skip-cell-clustering",
        action="store_true",
        help="Skip Scanpy QC/clustering; requires --celltype-key in table.obs",
    )
    parser.add_argument("--celltype-key", default="leiden", help="Obs column for Sainsc")
    parser.add_argument("--cell-type-col", default="cell_type", help="Obs column for cell-type tools")
    parser.add_argument(
        "--transcript-gene-source",
        default="feature_name",
        help="Transcript column copied to 'gene' when 'gene' is missing (Xenium / spatialdata-io)",
    )
    parser.add_argument("--binsize", type=float, default=5.0)
    parser.add_argument("--background-filter", type=float, default=0.2)
    parser.add_argument("--gaussian-kernel-key", type=float, default=2.0)
    parser.add_argument("--prob-threshold", type=float, default=0.5)
    parser.add_argument("--urna-square-size", type=float, default=50.0)
    parser.add_argument("--density-radius", type=float, default=10.0)
    parser.add_argument("--urna-source-max-dist", type=float, default=100.0)
    parser.add_argument("--urna-source-lmbda", type=float, default=0.1)
    parser.add_argument("--filter-min-logfc", type=float, default=7.0)
    parser.add_argument("--filter-max-p-val", type=float, default=0.05)
    parser.add_argument("--decontam-n-iter", type=int, default=5)
    parser.add_argument("--min-genes", type=int, default=10, help="sc.pp.filter_cells min_genes")
    parser.add_argument("--min-counts", type=int, default=30, help="sc.pp.filter_cells min_counts")
    parser.add_argument("--min-cells", type=int, default=3, help="sc.pp.filter_genes min_cells")
    parser.add_argument(
        "--spatial-genes",
        default="",
        help="Comma-separated genes for tp.pl.spatial_transcripts (empty = skip)",
    )
    parser.add_argument(
        "--diffusion-genes",
        default="",
        help="Comma-separated genes for tp.pl.gene_distribution_from_source (empty = skip)",
    )
    parser.add_argument(
        "--source-spatial-genes",
        default="",
        help="Comma-separated genes for tp.pl.spatial_transcripts_source (empty = skip)",
    )
    parser.add_argument(
        "--roi",
        default="",
        help="Optional ROI as x_min,x_max,y_min,y_max for spatial_transcripts_source",
    )
    parser.add_argument(
        "--umap-vmax",
        default="",
        help="Optional vmax for UMAP urna_source_score plot (e.g. p99.995)",
    )
    parser.add_argument(
        "--contact-vmax-combined",
        type=float,
        default=None,
        help="Optional vmax for combined cell contact heatmap",
    )
    parser.add_argument(
        "--contact-vmax-urna",
        type=float,
        default=None,
        help="Optional vmax for uRNA-specific cell contact heatmap",
    )
    return parser.parse_args()


def comma_separated_list(value: str) -> list[str]:
    if not value or not value.strip():
        return []
    return [part.strip() for part in value.split(",") if part.strip()]


def parse_roi(value: str) -> tuple[float, float, float, float] | None:
    if not value or not value.strip():
        return None
    parts = [float(x.strip()) for x in value.split(",")]
    if len(parts) != 4:
        raise ValueError("--roi must be four comma-separated numbers: x_min,x_max,y_min,y_max")
    return parts[0], parts[1], parts[2], parts[3]


def save_figure(result: Any, plots_dir: Path, stem: str) -> None:
    import matplotlib.pyplot as plt

    plots_dir.mkdir(parents=True, exist_ok=True)
    out_path = plots_dir / f"{stem}.png"

    fig = None
    if isinstance(result, tuple) and result and hasattr(result[0], "savefig"):
        fig = result[0]
    elif result is not None and hasattr(result, "savefig"):
        fig = result
    else:
        candidate = plt.gcf()
        if candidate.get_axes():
            fig = candidate

    if fig is None:
        return

    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)


def prepare_transcript_gene_column(sdata, gene_key: str, source_col: str) -> None:
    """Ensure points layer has `gene` (spatialdata-io Xenium uses feature_name)."""
    transcripts = sdata["transcripts"]
    columns = list(transcripts.columns) if hasattr(transcripts, "columns") else []

    if gene_key in columns:
        print(f"Transcript column '{gene_key}' already present.", file=sys.stderr)
    elif source_col in columns:
        print(f"Creating '{gene_key}' from '{source_col}'.", file=sys.stderr)
        transcripts[gene_key] = transcripts[source_col]
    else:
        raise ValueError(
            f"Transcripts need '{gene_key}' or '{source_col}'; found: {columns[:20]}"
        )

    transcripts[gene_key] = transcripts[gene_key].astype("category")

    # Sainsc needs 'overlaps_cell' (normally set by tp.pp.xenium_converter)
    if "overlaps_cell" not in columns:
        if "cell_id" not in columns:
            raise ValueError("Transcripts need 'overlaps_cell' or 'cell_id' to derive it")
        print("Creating 'overlaps_cell' from 'cell_id' != 'UNASSIGNED'.", file=sys.stderr)
        transcripts["overlaps_cell"] = transcripts["cell_id"].astype(str) != "UNASSIGNED"


def ensure_raw_layer(adata) -> None:
    """Decontamination uses adata.layers['raw'], not adata.raw alone."""
    if "raw" in adata.layers:
        return
    if adata.raw is not None:
        adata.layers["raw"] = adata.raw.X.copy()
    else:
        adata.layers["raw"] = adata.X.copy()
    print("Created layers['raw'] from expression matrix.", file=sys.stderr)


def ensure_umap(adata) -> None:
    import scanpy as sc

    if "X_umap" in adata.obsm:
        return
    print("Computing UMAP for plotting...", file=sys.stderr)
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata)
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=10)
    sc.tl.umap(adata)


def prepare_cell_table(
    sdata,
    celltype_key: str,
    cell_type_col: str,
    skip_cell_clustering: bool,
    need_umap: bool,
    min_genes: int = 10,
    min_counts: int = 30,
    min_cells: int = 3,
) -> None:
    """Tutorial §3: QC, normalize, cluster; leiden + cell_type + layers['raw']."""
    import scanpy as sc

    adata = sdata["table"]
    adata.raw = adata.copy()

    if skip_cell_clustering:
        if celltype_key not in adata.obs.columns:
            raise ValueError(
                f"--skip-cell-clustering set but obs has no '{celltype_key}' column"
            )
        ensure_raw_layer(adata)
        if need_umap:
            ensure_umap(adata)
    else:
        print(f"Running Scanpy preprocessing and Leiden ('{celltype_key}')...", file=sys.stderr)
        all_cells = adata.shape[0]
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_cells(adata, min_counts=min_counts)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        print(
            f"Cells retained after filtering: {adata.shape[0] / max(all_cells, 1):.2%}",
            file=sys.stderr,
        )
        adata.layers["raw"] = adata.X.copy()
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata, n_neighbors=10)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, key_added=celltype_key, resolution=0.6)
        sc.tl.rank_genes_groups(adata, groupby=celltype_key, use_raw=False)

    if cell_type_col not in adata.obs.columns:
        print(
            f"Column '{cell_type_col}' missing; using string labels from '{celltype_key}'.",
            file=sys.stderr,
        )
        adata.obs[cell_type_col] = adata.obs[celltype_key].astype(str)

    ensure_raw_layer(adata)


def export_tables(
    sdata,
    tables_dir: Path,
    spillover,
) -> None:
    tables_dir.mkdir(parents=True, exist_ok=True)
    if "xrna_metadata" in sdata.tables:
        var = sdata["xrna_metadata"].var
        var.to_csv(tables_dir / "xrna_metadata_var.csv")
        if "moran_I" in var.columns:
            var.sort_values(by="moran_I", ascending=False).head(50).to_csv(
                tables_dir / "top_spatially_variable_genes.csv"
            )
    if spillover is not None:
        spillover.to_csv(tables_dir / "intracellular_spillover.csv")
    if "source_score" in sdata.tables:
        score = sdata["source_score"]
        out = tables_dir / "source_score.csv"
        if hasattr(score, "obs"):
            score.obs.to_csv(out)
        elif hasattr(score, "to_pandas"):
            score.to_pandas().to_csv(out)
        elif hasattr(score, "compute"):
            score.compute().to_csv(out)


def run_analysis(args: argparse.Namespace) -> None:
    import scanpy as sc
    import spatialdata as sd
    import troutpy as tp

    if not args.skip_plots:
        import matplotlib

        matplotlib.use("Agg")

    spatial_genes = comma_separated_list(args.spatial_genes)
    diffusion_genes = comma_separated_list(args.diffusion_genes)
    source_spatial_genes = comma_separated_list(args.source_spatial_genes)
    roi = parse_roi(args.roi)

    print(f"Reading SpatialData from {args.input}", file=sys.stderr)
    sdata = sd.read_zarr(args.input)

    prepare_transcript_gene_column(sdata, "gene", args.transcript_gene_source)
    prepare_cell_table(
        sdata,
        args.celltype_key,
        args.cell_type_col,
        args.skip_cell_clustering,
        need_umap=not args.skip_plots,
        min_genes=args.min_genes,
        min_counts=args.min_counts,
        min_cells=args.min_cells,
    )

    print("Segmentation-free Sainsc...", file=sys.stderr)
    tp.pp.segmentation_free_sainsc(
        sdata,
        binsize=args.binsize,
        celltype_key=args.celltype_key,
        background_filter=args.background_filter,
        gaussian_kernel_key=args.gaussian_kernel_key,
    )

    if not args.skip_plots:
        save_figure(
            tp.pl.histogram(
                sdata,
                x="cosine_similarity",
                hue="overlaps_cell",
                palette="troutpy",
            ),
            args.plots_dir,
            "histogram_cosine_similarity",
        )

    print("Define uRNA (extracellular)...", file=sys.stderr)
    tp.pp.define_urna(
        sdata,
        layer="transcripts",
        method="sainsc",
        prob_threshold=args.prob_threshold,
    )

    if not args.skip_plots:
        save_figure(
            tp.pl.crosstab(
                sdata,
                yvar="extracellular",
                xvar="overlaps_cell",
                normalize=False,
                cmap="troutpy",
                kind="barh",
                stacked=True,
                figsize=(5, 2),
            ),
            args.plots_dir,
            "crosstab_extracellular_overlaps_cell",
        )

    print("Gene-centric uRNA metrics...", file=sys.stderr)
    tp.tl.quantify_overexpression(
        sdata,
        layer="transcripts",
        codeword_key="codeword_category",
        control_codewords=DEFAULT_CONTROL_CODEWORDS,
        gene_key="gene",
    )

    if not args.skip_plots:
        save_figure(
            tp.pl.metric_scatter(
                sdata,
                x="count",
                y="logfoldchange_over_noise",
                label_top_n_x=7,
                label_top_n_y=5,
                label_bottom_n_x=0,
                label_bottom_n_y=0,
                size=20,
            ),
            args.plots_dir,
            "metric_scatter_count_vs_logfc_noise",
        )

    tp.tl.extracellular_enrichment(sdata)

    if not args.skip_plots:
        save_figure(
            tp.pl.metric_scatter(
                sdata,
                x="extracellular_proportion",
                y="logfoldchange_over_noise",
                label_top_n_x=7,
                label_top_n_y=5,
                label_bottom_n_x=0,
                label_bottom_n_y=0,
                size=5,
                linewidth=0.1,
            ),
            args.plots_dir,
            "metric_scatter_extracellular_prop_vs_logfc",
        )

    print("Spatial variability and in/out correlation...", file=sys.stderr)
    tp.pp.aggregate_urna(sdata, gene_key="gene", square_size=args.urna_square_size)
    tp.tl.spatial_variability(sdata, gene_key="gene", n_neighbors=10)
    tp.tl.in_out_correlation(sdata, n_neighbors=20)

    if not args.skip_plots:
        for gene in spatial_genes:
            save_figure(
                tp.pl.spatial_transcripts(
                    sdata,
                    gene_list=[gene],
                    scatter_size=2,
                    boundary_linewidth=0.1,
                ),
                args.plots_dir,
                f"spatial_transcripts_{gene}",
            )
        save_figure(
            tp.pl.metric_scatter(
                sdata,
                x="moran_I",
                y="in_out_spearmanR",
                label_top_n_x=0,
                label_top_n_y=3,
                label_bottom_n_x=0,
                label_bottom_n_y=3,
                size=5,
                linewidth=0.1,
            ),
            args.plots_dir,
            "metric_scatter_moran_vs_in_out_spearman",
        )

    print("Filter uRNA genes...", file=sys.stderr)
    tp.pp.filter_urna(
        sdata,
        gene_key="gene",
        min_logfoldratio_over_noise=args.filter_min_logfc,
        max_p_val_noise=args.filter_max_p_val,
        control_probe=False,
    )

    print("uRNA source and diffusion...", file=sys.stderr)
    tp.tl.density_similarity(sdata, radius=args.density_radius)
    tp.tl.urna_source_score(
        sdata,
        cell_type_col=args.cell_type_col,
        max_dist=args.urna_source_max_dist,
        lmbda=args.urna_source_lmbda,
    )
    tp.tl.assess_diffusion(sdata, gene_key="gene", distance_key="distance_to_source")

    if not args.skip_plots:
        save_figure(
            tp.pl.diffusion_results(sdata, label_top_n_y=10, y_logscale=True),
            args.plots_dir,
            "diffusion_results",
        )
        if diffusion_genes:
            save_figure(
                tp.pl.gene_distribution_from_source(
                    sdata,
                    diffusion_genes,
                    distance_key="distance_to_source",
                ),
                args.plots_dir,
                "gene_distribution_from_source",
            )
        save_figure(
            tp.pl.source_score_by_celltype(sdata, figsize=(6, 6), cmap="terrain"),
            args.plots_dir,
            "source_score_by_celltype",
        )
        umap_kwargs: dict[str, Any] = {
            "color": ["urna_source_score", args.cell_type_col],
            "cmap": "terrain",
            "title": "uRNA source score",
            "show": False,
        }
        if args.umap_vmax.strip():
            umap_kwargs["vmax"] = args.umap_vmax.strip()
        save_figure(
            sc.pl.umap(sdata["table"], **umap_kwargs),
            args.plots_dir,
            "umap_urna_source_score",
        )
        if source_spatial_genes:
            source_kwargs: dict[str, Any] = {
                "cell_color_key": args.cell_type_col,
                "extra_color_key": args.cell_type_col,
                "gene_list": source_spatial_genes,
                "point_size": 3,
                "show_edges": True,
            }
            if roi is not None:
                source_kwargs["use_roi"] = True
                source_kwargs["roi"] = roi
            save_figure(
                tp.pl.spatial_transcripts_source(sdata, **source_kwargs),
                args.plots_dir,
                "spatial_transcripts_source",
            )

    print("Target scores and cell contacts...", file=sys.stderr)
    tp.tl.compute_target_score(
        sdata,
        layer="transcripts",
        lambda_decay=0.1,
        k_neighbors=min(50, sdata["table"].n_obs),  # tiny test data has <50 cells
        copy=False,
        celltype_key=args.cell_type_col,
        gene_key="gene",
    )
    tp.tl.cell_contacts_with_urna_sources(
        sdata,
        spatial_key="spatial",
        cell_type_key=args.cell_type_col,
        distance=20,
    )

    if not args.skip_plots:
        combined_kwargs: dict[str, Any] = {
            "celltype_key": args.cell_type_col,
            "figsize": (7, 7),
            "dendrogram_ratio": 0.1,
            "cmap": "YlGnBu",
            "key": "cell_contact_combined",
        }
        if args.contact_vmax_combined is not None:
            combined_kwargs["vmax"] = args.contact_vmax_combined
        save_figure(
            tp.pl.cell_type_contacts(sdata, **combined_kwargs),
            args.plots_dir,
            "cell_type_contacts_combined",
        )
        urna_kwargs: dict[str, Any] = {
            "celltype_key": args.cell_type_col,
            "figsize": (7, 7),
            "dendrogram_ratio": 0.1,
            "cmap": "YlGnBu",
            "key": "cell_contact_urna_specific",
        }
        if args.contact_vmax_urna is not None:
            urna_kwargs["vmax"] = args.contact_vmax_urna
        save_figure(
            tp.pl.cell_type_contacts(sdata, **urna_kwargs),
            args.plots_dir,
            "cell_type_contacts_urna_specific",
        )

    print("Decontamination...", file=sys.stderr)
    spillover = tp.tl.flag_intracellular_spillover(
        sdata,
        prob_threshold=args.prob_threshold,
        cell_type_col=args.celltype_key,
    )
    tp.tl.decontaminate_cell_expression(
        sdata,
        spillover,
        layer_key="raw",
        output_layer="decontaminated",
    )
    tp.tl.subtract_local_urna_background(
        sdata,
        radius=100,
        inner_radius=10,
        raw_layer="raw",
        output_layer="decontaminated_density",
    )

    adata = sdata["table"]
    for layer in ["raw", "decontaminated", "decontaminated_density"]:
        if layer in adata.layers:
            print(f"{layer:>24}: {adata.layers[layer].sum():,.0f} counts", file=sys.stderr)

    tp.tl.iterative_soft_decontamination(
        sdata,
        n_iter=args.decontam_n_iter,
        raw_layer="raw",
        extracellular_score_kwargs={"cell_type_col": args.cell_type_col},
    )

    export_tables(sdata, args.tables_dir, spillover)

    args.output_zarr.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing SpatialData to {args.output_zarr}", file=sys.stderr)
    sdata.write(args.output_zarr)


def main() -> None:
    args = parse_args()
    if not args.input.is_dir():
        print(f"Input Zarr not found or not a directory: {args.input}", file=sys.stderr)
        sys.exit(1)
    run_analysis(args)


if __name__ == "__main__":
    main()
