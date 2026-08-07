#!/usr/bin/env python3

"""
Transcript QC Processing Module
Performs all analysis from notebooks/1_qc_molecule.ipynb and generates figures and metrics.
Authors: Malwina Prater, mprater@altoslabs.com; Dongze He, dhe@altoslabs.com; Felix Krueger, fkrueger@altoslabs.com
"""

import argparse
import json
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import pyarrow.parquet as pq


def calculate_noise_bound(n_molecules_non_gene_prefix, quant: float = 0.99):
    """Calculate the noise bounds based on the non-gene molecules."""
    from scipy.stats import median_abs_deviation, norm

    if n_molecules_non_gene_prefix.empty:
        return 0, 0

    quant_val = norm.ppf(quant)
    n_mols_log = np.log10(n_molecules_non_gene_prefix.values)
    std = median_abs_deviation(n_mols_log, scale="normal")
    noise_lb = np.mean(n_mols_log) - quant_val * std
    noise_ub = np.mean(n_mols_log) + quant_val * std
    return 10**noise_lb, 10**noise_ub


def estimate_min_mols_per_cell(n_mols_per_cell, min_value: int = 10):
    n_mols_per_cell = np.log10(np.asarray(n_mols_per_cell) + 1)
    nm_hist = np.histogram(n_mols_per_cell, bins=100)
    mode = nm_hist[1][nm_hist[0].argmax()]
    ci = np.quantile(n_mols_per_cell[n_mols_per_cell > mode], 0.99) - mode
    return max(min_value, int(round(10 ** (mode - ci))))


# Set plotting style
sns.set_theme(style="whitegrid")
plt.rcParams["figure.dpi"] = 300


def write_versions(outdir: Path, task_process: str) -> None:
    """Write versions.yml into the output dir for the transcript QC report to display.

    This is separate from the module's topic-channel version reporting: it records
    the analysis package versions for the rendered HTML report. A missing package
    raises (importlib.metadata.version), so a broken environment surfaces loudly.
    """
    from importlib.metadata import version
    import platform

    packages = [
        "matplotlib",
        "numpy",
        "pandas",
        "scipy",
        "seaborn",
        "scanpy",
        "pyarrow",
    ]
    lines = [f'"{task_process}":']
    for pkg in packages:
        lines.append(f"    {pkg}: {version(pkg)}")
    lines.append(f"    python: {platform.python_version()}")
    (outdir / "versions.yml").write_text("\n".join(lines) + "\n")


def read_random_parquet_row_groups(parquet_file, num_row_groups=4, random_seed=42):
    """
    Read a random subset of row groups from a parquet file. Each row group is a set of rows that are contiguous in the file. For 10x transcripts.parquet file, each row group has about 262,000 rows.
    """
    np.random.seed(random_seed)
    selected_row_groups = np.random.choice(
        parquet_file.metadata.num_row_groups, size=num_row_groups, replace=False
    )
    return parquet_file.read_row_groups(selected_row_groups).to_pandas()


def main():
    parser = argparse.ArgumentParser(description="Transcript QC Processing")
    parser.add_argument(
        "--xenium-bundle-dir", required=True, help="Path to Xenium bundle directory"
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--non-gene-prefix",
        default="NegControlProbe",
        help="Prefix for non-gene features",
    )
    parser.add_argument(
        "--stain-names", help="Stain names (unused but kept for compatibility)"
    )
    parser.add_argument(
        "--task-process", default="TRANSCRIPT_QC", help="Task process name"
    )
    parser.add_argument(
        "--num-row-groups",
        type=int,
        default=None,
        help="Number of row groups to process",
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Number of threads (for compatibility)"
    )

    args = parser.parse_args()

    # Validate parameters
    if args.xenium_bundle_dir is None or not os.path.exists(args.xenium_bundle_dir):
        raise FileNotFoundError(
            f'The given XENIUM_BUNDLE_DIR, "{args.xenium_bundle_dir}" doesn\'t exist'
        )

    XENIUM_BUNDLE_DIR = Path(args.xenium_bundle_dir)
    # Create output directories
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_fig_dir = outdir / "figures"
    output_fig_dir.mkdir(parents=True, exist_ok=True)
    figures_source_dir = outdir / "figures_source"
    figures_source_dir.mkdir(parents=True, exist_ok=True)
    output_metrics_path = outdir / "transcript_qc_metrics.json"

    NUM_ROW_GROUPS = args.num_row_groups

    transcripts_parquet_path = XENIUM_BUNDLE_DIR / "transcripts.parquet"
    morphology_focus_dir = XENIUM_BUNDLE_DIR / "morphology_focus"
    cells_parquet_path = XENIUM_BUNDLE_DIR / "cells.parquet"
    cell_feature_matrix_h5_path = XENIUM_BUNDLE_DIR / "cell_feature_matrix.h5"

    # EXACT CODE FROM ORIGINAL NOTEBOOK - check if all required files are present
    required_files = [
        transcripts_parquet_path,
        morphology_focus_dir,
        cell_feature_matrix_h5_path,
        cells_parquet_path,
    ]
    for file in required_files:
        if not os.path.exists(file):
            print(f"Required file not found: {file}")
            sys.exit(1)

    print("=== Transcript QC Processing ===")
    print(f"Input directory: {XENIUM_BUNDLE_DIR}")
    print(f"Output directory: {outdir}")

    # EXACT CODE FROM ORIGINAL NOTEBOOK - check if all required columns are present in the transcripts parquet file
    transcripts_parquet = pq.ParquetFile(transcripts_parquet_path)
    transcripts_parquet_columns = transcripts_parquet.schema.names
    num_molecules = transcripts_parquet.metadata.num_rows
    print("Total number of molecules: {:,}".format(num_molecules))

    # EXACT CODE FROM ORIGINAL NOTEBOOK - required columns
    required_columns = [
        "cell_id",
        "qv",
        "fov_name",
        "codeword_category",
        "is_gene",
        "feature_name",
    ]
    missing_columns = [
        col for col in required_columns if col not in transcripts_parquet_columns
    ]
    if missing_columns:
        print(f"Missing required columns in transcripts.parquet: {missing_columns}")
        sys.exit(1)

    # EXACT CODE FROM ORIGINAL NOTEBOOK - load data
    if NUM_ROW_GROUPS is None:
        df_spatial = pd.read_parquet(
            transcripts_parquet_path,
            columns=required_columns,
        )
    elif NUM_ROW_GROUPS > transcripts_parquet.metadata.num_row_groups:
        print(
            f"NUM_ROW_GROUPS ({NUM_ROW_GROUPS}) is greater than the number of molecules in the parquet file ({transcripts_parquet.metadata.num_row_groups}). Read the entire file."
        )
        NUM_ROW_GROUPS = None
        df_spatial = pd.read_parquet(transcripts_parquet_path, columns=required_columns)
    else:
        df_spatial = read_random_parquet_row_groups(transcripts_parquet, NUM_ROW_GROUPS)

    num_selected_molecules = df_spatial.shape[0]

    if num_selected_molecules != num_molecules:
        print(
            f"Number of random molecules selected for analysis: {num_selected_molecules:,} (out of {num_molecules:,})"
        )

    codeword_category_counts = df_spatial["codeword_category"].value_counts()

    print(f"Features: {df_spatial.feature_name.unique().size:,}")
    print("\nFeature categories:")
    for cc in codeword_category_counts.sort_values(ascending=False).keys():
        count = codeword_category_counts[cc]
        percentage = count / num_selected_molecules * 100
        print(f"  {cc:<26} - {count:>12,} molecules ({percentage:>6.3f}%)")

    # EXACT CODE FROM ORIGINAL NOTEBOOK - quality distribution plots
    print("\nGenerating quality distribution plots...")

    # Filter data for quality analysis
    num_gene_molecules = df_spatial["is_gene"].sum()
    df_spatial_nongene = df_spatial.query("is_gene == False").sample(
        min(1000000, num_selected_molecules - num_gene_molecules), random_state=42
    )
    df_spatial_gene = df_spatial.query("is_gene == True").sample(
        min(1000000, num_gene_molecules), random_state=42
    )

    # define hue order: genes first, then non-genes
    codeword_categories_order = codeword_category_counts.index
    codeword_categories_order = (
        codeword_categories_order[codeword_categories_order.str.endswith("gene")]
        .sort_values(ascending=False)
        .tolist()
        + codeword_categories_order[~codeword_categories_order.str.endswith("gene")]
        .sort_values()
        .tolist()
    )

    df_spatial_quality = pd.concat([df_spatial_nongene, df_spatial_gene])

    print(
        f"Using {len(df_spatial_nongene):,} non-gene molecules and {len(df_spatial_gene):,} gene molecules for quality values (qv) distribution"
    )

    # Print total number of rows and rows per category in df_spatial_quality
    print("\nmolecule per category:")
    print(df_spatial_quality["codeword_category"].value_counts())

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Figure 1: Quality distribution density plot with overlapping curves
    fig = plt.figure(figsize=(15, 8))

    # Create density plots with all curves in the same plot
    for i, category in enumerate(codeword_categories_order):
        # Filter data for this category
        category_data = df_spatial_quality[
            df_spatial_quality["codeword_category"] == category
        ]["qv"]

        # Create density plot for this category
        sns.kdeplot(x=category_data, fill=True, alpha=0.5, linewidth=2, label=category)

    plt.xlabel("Quality Value (qv)", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.title(
        "Distribution of Transcript Quality by Codeword Category", fontsize=14, pad=20
    )

    # Add vertical line at qv=20
    plt.axvline(
        x=20, color="darkred", linestyle="--", linewidth=2, label="QV threshold = 20"
    )

    # Add legend
    plt.legend(title="Codeword Category", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(
        output_fig_dir / "quality_distribution_density.pdf",
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        output_fig_dir / "quality_distribution_density.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Save data for quality distribution density plot
    df_spatial_quality.to_csv(
        figures_source_dir / "quality_distribution_density.csv", index=False
    )

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Figure 2: Quality distribution by codeword category (violin plot)
    fig = plt.figure(figsize=(15, 8))
    sns.violinplot(
        data=df_spatial_quality,
        x="codeword_category",
        y="qv",
        hue="codeword_category",
        split=False,
        inner="box",
        palette="husl",
        density_norm="width",
        legend=False,
        order=codeword_categories_order,
    )
    plt.xlabel("Codeword Category", fontsize=12)
    plt.ylabel("Quality Value (qv)", fontsize=12)
    plt.title(
        "Distribution of Transcript Quality by Codeword Category", fontsize=14, pad=20
    )
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    # add a horizontal line at qv=20
    plt.axhline(y=20, color="grey", linestyle="--")
    plt.savefig(
        output_fig_dir / "quality_distributions_comprehensive.pdf",
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        output_fig_dir / "quality_distributions_comprehensive.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Save data for quality distributions comprehensive
    df_spatial_quality.to_csv(
        figures_source_dir / "quality_distributions_comprehensive.csv", index=False
    )

    # Save df_spatial_quality to CSV
    output_file = outdir / "df_spatial_quality.csv"
    df_spatial_quality.to_csv(output_file, index=False)
    print(f"Saved df_spatial_quality to: {output_file}")

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Figure 3: Quality distribution by Field of View
    fig = plt.figure(figsize=(15, 8))
    sns.violinplot(
        data=df_spatial_gene,
        x="fov_name",
        y="qv",
        hue="codeword_category",
        split=False,
        inner="box",
        palette="husl",
        density_norm="width",
        legend=False,
    )
    plt.xlabel("Field of View", fontsize=12)
    plt.ylabel("Quality Value (qv)", fontsize=12)
    plt.title(
        "Distribution of Transcript Quality by Field of View", fontsize=14, pad=20
    )
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    # add a horizontal line at qv=20
    plt.axhline(y=20, color="grey", linestyle="--")
    plt.savefig(output_fig_dir / "quality_by_fov.pdf", dpi=300, bbox_inches="tight")
    plt.savefig(output_fig_dir / "quality_by_fov.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Save data for quality by fov
    df_spatial_gene.to_csv(figures_source_dir / "quality_by_fov.csv", index=False)

    del df_spatial_quality, df_spatial_gene

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Calculate noise threshold
    _, n_mols_threshold = calculate_noise_bound(
        df_spatial_nongene["feature_name"][
            df_spatial_nongene["feature_name"].str.startswith("NegControl")
        ].value_counts()
    )
    n_mols_threshold = n_mols_threshold * (num_molecules / num_selected_molecules)
    print(
        f"Noise threshold for genes' molecule count: {n_mols_threshold:.0f} molecules"
    )

    # EXACT CODE FROM ORIGINAL NOTEBOOK - group by feature_name
    n_mols_per_gene_df = df_spatial.groupby("feature_name").agg(
        # count of molecules per gene
        n_molecules=("feature_name", "count"),
        is_gene=("is_gene", "first"),
    )
    # rescale the number of molecules to account for subsampling
    n_mols_per_gene_df["n_molecules"] = n_mols_per_gene_df["n_molecules"] * (
        num_molecules / num_selected_molecules
    )

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Figure 4: Distribution of molecules per feature
    fig = plt.figure(figsize=(8, 4))
    sns.histplot(
        n_mols_per_gene_df,
        x="n_molecules",
        multiple="layer",  # Overlap instead of stack
        hue="is_gene",
        log_scale=True,
        bins=50,
        ax=plt.gca(),
        element="step",  # Outlined bars
        fill=False,  # No fill, just lines
    )
    plt.xlabel("Num. molecules")
    plt.ylabel("Num. features")
    plt.axvline(x=n_mols_threshold, color="grey", linestyle="--")
    plt.title("Distribution of molecules per feature", fontsize=14, pad=20)
    plt.tight_layout()
    plt.savefig(
        f"{output_fig_dir}/num_transcripts_per_feature.pdf",
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        f"{output_fig_dir}/num_transcripts_per_feature.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Save data for molecules per feature
    df_molecules_per_feature = n_mols_per_gene_df.copy()
    df_molecules_per_feature["n_mols_threshold"] = n_mols_threshold
    df_molecules_per_feature.to_csv(
        figures_source_dir / "num_transcripts_per_feature.csv", index=True
    )

    retained_genes = n_mols_per_gene_df.query(
        "n_molecules > @n_mols_threshold and is_gene == True"
    )
    print(
        f"Number of genes with a total molecule count higher than the threshold : {len(retained_genes):,}"
    )
    retained_genes.to_csv(outdir / "retained_genes.csv", index=True)

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Cell size distribution
    print("\nGenerating cell statistics plots...")

    # Check available columns and use appropriate cell area column
    available_columns = pq.ParquetFile(cells_parquet_path).schema.names
    cell_area_column = "cell_area" if "cell_area" in available_columns else "volume"
    cells_parquet = pd.read_parquet(
        cells_parquet_path, columns=["cell_id", cell_area_column]
    )
    # Rename the column for consistency
    cells_parquet.rename(columns={cell_area_column: "cell_size"}, inplace=True)

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Figure 5: Cell size distribution
    fig = plt.figure(figsize=(8, 6))
    sns.kdeplot(
        data=cells_parquet, x="cell_size", fill=True, color="skyblue", alpha=0.5
    )
    plt.xlabel("Genes per cell")
    plt.ylabel("Density")
    plt.title("Cell size distribution")

    # Add vertical lines for mean and median
    plt.axvline(
        cells_parquet["cell_size"].mean(), color="red", linestyle="--", label="Mean"
    )
    plt.axvline(
        cells_parquet["cell_size"].median(),
        color="green",
        linestyle="--",
        label="Median",
    )
    plt.legend()

    # Save the plot
    plt.savefig(
        os.path.join(outdir, "figures", "genes_per_cell_distribution.pdf"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        os.path.join(outdir, "figures", "genes_per_cell_distribution.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Save data for cell size distribution
    df_cell_size = cells_parquet[["cell_id", "cell_size"]].copy()
    df_cell_size["mean_cell_size"] = cells_parquet["cell_size"].mean()
    df_cell_size["median_cell_size"] = cells_parquet["cell_size"].median()
    df_cell_size.to_csv(
        figures_source_dir / "genes_per_cell_distribution.csv", index=False
    )
    del cells_parquet

    # Figure 6: Nucleus RNA fraction per cell
    cells_parquet = pd.read_parquet(
        cells_parquet_path,
        columns=["nucleus_count", "total_counts"],
        filters=[("total_counts", ">", 0)],
    )
    nucleus_count_fraction = cells_parquet["nucleus_count"] / (
        cells_parquet["total_counts"] + 1
    )

    # Create density plot if requested
    fig = plt.figure(figsize=(8, 6))
    sns.kdeplot(x=nucleus_count_fraction, fill=True, color="skyblue", alpha=0.5)
    plt.xlabel("Nucleus molecule fraction per cell")
    plt.ylabel("Density")
    # if nucleus_count_fraction is all 0, set the title to "nucleus molecule fraction distribution"
    plt.title(
        "no nucleus molecule detected; Empty plot shown"
        if nucleus_count_fraction.all() == 0
        else "nucleus molecule fraction distribution"
    )

    # Save the plot
    plt.savefig(
        os.path.join(
            outdir, "figures", "nucleus_transcript_fraction_per_cell_distribution.pdf"
        ),
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        os.path.join(
            outdir, "figures", "nucleus_transcript_fraction_per_cell_distribution.png"
        ),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Save data for nucleus molecule fraction
    df_nucleus_fraction = pd.DataFrame(
        {
            "nucleus_count": cells_parquet["nucleus_count"],
            "total_counts": cells_parquet["total_counts"],
            "nucleus_fraction": nucleus_count_fraction,
        }
    )
    df_nucleus_fraction.to_csv(
        figures_source_dir / "nucleus_transcript_fraction_per_cell_distribution.csv",
        index=False,
    )
    del cells_parquet

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Figure 7: Nucleus-to-cell area fraction
    cells_parquet = pd.read_parquet(
        cells_parquet_path, columns=["nucleus_area", "cell_area"]
    )
    nucleus_size_fraction = cells_parquet["nucleus_area"] / (
        cells_parquet["cell_area"] + 1
    )

    # Create density plot if requested
    fig = plt.figure(figsize=(8, 6))
    sns.kdeplot(x=nucleus_size_fraction, fill=True, color="skyblue", alpha=0.5)
    plt.xlabel("Nucleus to cell fraction")
    plt.ylabel("Density")
    plt.title(
        "nucleus to cell fraction distribution"
        if nucleus_size_fraction.all() != 0
        else "no nucleus molecule detected; Empty plot shown"
    )

    # Save the plot
    plt.savefig(
        os.path.join(
            outdir, "figures", "nucleus_to_cell_size_fraction_per_cell_distribution.pdf"
        ),
        dpi=300,
        bbox_inches="tight",
    )
    plt.savefig(
        os.path.join(
            outdir, "figures", "nucleus_to_cell_size_fraction_per_cell_distribution.png"
        ),
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)

    # Save data for nucleus to cell size fraction
    df_nucleus_size_fraction = pd.DataFrame(
        {
            "nucleus_area": cells_parquet["nucleus_area"],
            "cell_area": cells_parquet["cell_area"],
            "nucleus_size_fraction": nucleus_size_fraction,
        }
    )
    df_nucleus_size_fraction.to_csv(
        figures_source_dir / "nucleus_to_cell_size_fraction_per_cell_distribution.csv",
        index=False,
    )
    del cells_parquet

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Load cell feature matrix
    ad = sc.read_10x_h5(cell_feature_matrix_h5_path)
    # filter for retained genes
    ad = ad[:, ad.var_names.isin(retained_genes.index)]
    print(f"AnnData object with n_obs × n_vars = {ad.shape[0]} × {ad.shape[1]}")

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Figure 8: Distribution of molecules per cell
    n_mols_per_cell = ad.X.sum(axis=1).A1
    n_mols_threshold_cell = estimate_min_mols_per_cell(n_mols_per_cell)

    fig = plt.figure(figsize=(8, 4))
    sns.histplot(n_mols_per_cell, log_scale=True, bins=50, ax=plt.gca())
    plt.xlabel("Num. molecules")
    plt.ylabel("Num. cells")
    plt.axvline(x=n_mols_threshold_cell, color="grey", linestyle="--")
    plt.title("Distribution of molecules per Cell", fontsize=14, pad=20)
    plt.tight_layout()
    plt.savefig(
        f"{output_fig_dir}/num_transcripts_per_cell.pdf", dpi=300, bbox_inches="tight"
    )
    plt.savefig(
        f"{output_fig_dir}/num_transcripts_per_cell.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)
    print(f"Threshold for molecules per cell: {n_mols_threshold_cell}")

    # Save data for molecules per cell
    df_molecules_per_cell = pd.DataFrame(
        {
            "n_molecules_per_cell": n_mols_per_cell,
            "n_mols_threshold_cell": n_mols_threshold_cell,
        }
    )
    df_molecules_per_cell.to_csv(
        figures_source_dir / "num_transcripts_per_cell.csv", index=False
    )

    # Convert numpy array to pandas DataFrame
    n_mols_per_cell_df = pd.DataFrame(n_mols_per_cell, columns=["num_of_molecules"])
    # Save to CSV using os.path.join for path handling
    output_file = os.path.join(outdir, "num_transcripts_per_cell.csv")
    n_mols_per_cell_df.to_csv(output_file, index=False)
    print(f"Saved transcript distribution to: {output_file}")

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Figure 9: Distribution of genes per cell
    n_genes_per_cell = (ad.X != 0).sum(axis=1).A1
    n_genes_threshold = estimate_min_mols_per_cell(n_genes_per_cell)

    fig = plt.figure(figsize=(8, 4))
    sns.histplot(n_genes_per_cell, log_scale=True, bins=50, ax=plt.gca())
    plt.xlabel("Num. genes")
    plt.ylabel("Num. cells")
    plt.axvline(x=n_genes_threshold, color="grey", linestyle="--")
    plt.title("Distribution of number of detected genes per Cell", fontsize=14, pad=20)
    plt.tight_layout()
    plt.savefig(
        f"{output_fig_dir}/num_genes_per_cell.pdf", dpi=300, bbox_inches="tight"
    )
    plt.savefig(
        f"{output_fig_dir}/num_genes_per_cell.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)

    # Save data for genes per cell
    df_genes_per_cell = pd.DataFrame(
        {"n_genes_per_cell": n_genes_per_cell, "n_genes_threshold": n_genes_threshold}
    )
    df_genes_per_cell.to_csv(figures_source_dir / "num_genes_per_cell.csv", index=False)

    # Convert numpy array to pandas DataFrame
    n_genes_per_cell_df = pd.DataFrame(n_genes_per_cell, columns=["num_of_genes"])
    # Save to CSV using os.path.join for path handling
    output_file = os.path.join(outdir, "num_genes_per_cell.csv")
    n_genes_per_cell_df.to_csv(output_file, index=False)
    print(f"Saved gene distribution to: {output_file}")

    # EXACT CODE FROM ORIGINAL NOTEBOOK - Save metrics
    metrics = {
        "total_transcripts": int(num_molecules),
        "selected_transcripts": int(num_selected_molecules),
        "total_features": int(df_spatial.feature_name.nunique()),
        "codeword_category_counts": {
            str(k): int(v) for k, v in codeword_category_counts.items()
        },
        "neg_control_quantile": int(n_mols_threshold),
        "min_transcripts_per_cell": int(n_mols_threshold_cell),
        "min_genes_per_cell": int(n_genes_threshold),
        "retained_genes_count": len(retained_genes),
        "total_cells": int(ad.shape[0]),
        "analyzed_genes": int(ad.shape[1]),
    }

    # Save metrics
    with open(output_metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    # Save versions.yml — read and displayed by the transcript QC report
    write_versions(outdir, args.task_process)

    print("\n=== Processing Complete ===")
    print(f"Generated {len(list(output_fig_dir.glob('*.pdf')))} figures")
    print(f"Generated {len(list(figures_source_dir.glob('*.csv')))} CSV data files")
    print(f"Saved metrics to: {output_metrics_path}")
    print(f"Output directory: {outdir}")
    print(f"Figures directory: {output_fig_dir}")
    print(f"Source data directory: {figures_source_dir}")


if __name__ == "__main__":
    main()
