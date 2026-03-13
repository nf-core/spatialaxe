#!/usr/bin/env python3

import json
from pathlib import Path
import numpy as np
import pandas as pd
import tifffile

# Xenium full-resolution image: 1 pixel = 0.2125 µm (10x Genomics spec).
# Transcript x_location / y_location are in microns.
# To overlay transcripts on the full-res image: pixel = micron / pixel_size.
XENIUM_DEFAULT_PIXEL_SIZE_UM = 0.2125


def _pick_column(df: pd.DataFrame,
                 candidates: list[str],
                 required: bool = True):
    for name in candidates:
        if name in df.columns:
            return name
    if required:
        raise ValueError(f"Could not find any of the required columns: {candidates}")
    return None


def _read_pixel_size(experiment_xenium_path: str) -> float:
    """Read pixel_size (µm/px) from experiment.xenium; fall back to 0.2125."""
    try:
        with open(experiment_xenium_path) as fh:
            meta = json.load(fh)
        return float(meta.get("pixel_size", XENIUM_DEFAULT_PIXEL_SIZE_UM))
    except Exception:
        return XENIUM_DEFAULT_PIXEL_SIZE_UM


def convert_xenium_to_scs(parquet_path: str,
                          output_tsv: str,
                          output_bgi_tsv: str,
                          morphology_image_path: str,
                          output_morph2d_tif: str,
                          metrics_tsv: str,
                          experiment_xenium_path: str = "",
                          bin_size: float = 1.0):
    """
    Convert Xenium transcripts to SCS/BGI format with correct pixel-space coordinates.

    Xenium x_location / y_location are in microns.
    The morphology image full resolution is 0.2125 µm/px (Xenium spec).
    Coordinates are converted to pixels: pixel = micron / pixel_size.
    The morphology image is cropped to the pixel ROI covered by the transcripts.
    """
    pixel_size = _read_pixel_size(experiment_xenium_path) if experiment_xenium_path else XENIUM_DEFAULT_PIXEL_SIZE_UM

    transcripts = pd.read_parquet(parquet_path, engine="pyarrow")

    gene_col = _pick_column(transcripts, ["feature_name", "gene", "gene_id", "geneID"])
    x_col    = _pick_column(transcripts, ["x_location", "x", "x_global_px", "x_centroid"])
    y_col    = _pick_column(transcripts, ["y_location", "y", "y_global_px", "y_centroid"])
    count_col = _pick_column(transcripts, ["counts", "count", "n_counts"], required=False)

    table = transcripts[[gene_col, x_col, y_col]].copy()
    table = table.dropna(subset=[gene_col, x_col, y_col])

    # Convert micron coordinates → full-resolution pixel coordinates.
    # Xenium: x_location is along image width (columns), y_location along height (rows).
    table["row_px"]    = (table[y_col].astype(float) / pixel_size).round().astype(int)
    table["column_px"] = (table[x_col].astype(float) / pixel_size).round().astype(int)

    # Optionally merge into user-specified bins (bin_size in pixels, default 1 = no binning).
    if bin_size > 1:
        table["row"]    = (table["row_px"]    / bin_size).astype(int)
        table["column"] = (table["column_px"] / bin_size).astype(int)
    else:
        # Zero-base pixel coordinates so the BGI file starts at (0, 0).
        r0 = table["row_px"].min()
        c0 = table["column_px"].min()
        table["row"]    = table["row_px"]    - r0
        table["column"] = table["column_px"] - c0

    if count_col is None:
        table["counts"] = 1
    else:
        table["counts"] = transcripts.loc[table.index, count_col].fillna(1).astype(int)

    table = table.rename(columns={gene_col: "geneID"})[["geneID", "row", "column", "counts"]]
    table = table.groupby(["geneID", "row", "column"], as_index=False)["counts"].sum()

    out_tsv = Path(output_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_tsv, sep="\t", index=False)

    # ── Morphology image ────────────────────────────────────────────────────────
    # Load and collapse to 2D (max projection across z/channels).
    image = tifffile.imread(morphology_image_path)
    image = np.squeeze(np.asarray(image))
    if image.ndim == 2:
        image2d = image
    elif image.ndim >= 3:
        h, w = image.shape[-2], image.shape[-1]
        image2d = image.reshape((-1, h, w)).max(axis=0)
    else:
        raise ValueError(f"Unsupported morphology image shape: {image.shape}")

    # Crop to the pixel ROI covered by transcripts.
    # Derive absolute pixel bounds directly from physical coords in the parquet.
    r_min_abs = int(round(float(transcripts[y_col].min()) / pixel_size))
    r_max_abs = int(round(float(transcripts[y_col].max()) / pixel_size))
    c_min_abs = int(round(float(transcripts[x_col].min()) / pixel_size))
    c_max_abs = int(round(float(transcripts[x_col].max()) / pixel_size))

    # Clamp to image bounds.
    H, W = image2d.shape
    r_min_abs = max(0, r_min_abs)
    r_max_abs = min(H - 1, r_max_abs)
    c_min_abs = max(0, c_min_abs)
    c_max_abs = min(W - 1, c_max_abs)

    cropped = image2d[r_min_abs:r_max_abs + 1, c_min_abs:c_max_abs + 1]

    out_morph2d = Path(output_morph2d_tif)
    out_morph2d.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(out_morph2d, cropped)

    # ── BGI file (SCS/spateo format) ────────────────────────────────────────────
    # spateo read_bgi_agg: x → AnnData dim-0 (height/rows), y → dim-1 (width/cols).
    # Our table["row"] = height direction, table["column"] = width direction.
    bgi = pd.DataFrame({
        "geneID":    table["geneID"],
        "x":         table["row"].astype(int),
        "y":         table["column"].astype(int),
        "MIDCounts": table["counts"].astype(int),
    })

    out_bgi_tsv = Path(output_bgi_tsv)
    out_bgi_tsv.parent.mkdir(parents=True, exist_ok=True)
    bgi.to_csv(out_bgi_tsv, sep="\t", index=False)

    metrics = {
        "n_rows":         int(len(table)),
        "n_unique_genes": int(table["geneID"].nunique()),
        "row_min":        int(table["row"].min())    if len(table) else 0,
        "row_max":        int(table["row"].max())    if len(table) else 0,
        "column_min":     int(table["column"].min()) if len(table) else 0,
        "column_max":     int(table["column"].max()) if len(table) else 0,
        "pixel_size_um":  float(pixel_size),
        "bin_size":       float(bin_size),
        "morph2d_H":      int(cropped.shape[0]),
        "morph2d_W":      int(cropped.shape[1]),
    }

    pd.DataFrame(
        {"metric": list(metrics.keys()), "value": list(metrics.values())}
    ).to_csv(metrics_tsv, sep="\t", index=False)


if __name__ == "__main__":
    transcripts_parquet: str    = "${transcripts_parquet}"
    morphology_image: str       = "${morphology_image}"
    experiment_xenium: str      = "${experiment_xenium}"
    prefix: str                 = "${prefix}"
    bin_size: float             = float("${task.ext.bin_size ?: 1.0}")

    output_tsv        = f"{prefix}/scs_input.tsv"
    output_bgi_tsv    = f"{prefix}/scs_input_bgi.tsv"
    output_morph2d_tif = f"{prefix}/morph2d.tif"
    metrics_tsv       = f"{prefix}/xenium2scs_metrics.tsv"

    convert_xenium_to_scs(
        parquet_path=transcripts_parquet,
        output_tsv=output_tsv,
        output_bgi_tsv=output_bgi_tsv,
        morphology_image_path=morphology_image,
        output_morph2d_tif=output_morph2d_tif,
        metrics_tsv=metrics_tsv,
        experiment_xenium_path=experiment_xenium,
        bin_size=bin_size,
    )

    with open("versions.yml", "w", encoding="utf-8") as fobj:
        fobj.write('"${task.process}":\\n')
        fobj.write('xenium2scs: "1.0.0"\\n')
