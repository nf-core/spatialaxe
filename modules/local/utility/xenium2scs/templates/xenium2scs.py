#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import tifffile

# Xenium full-resolution image: 1 pixel = 0.2125 µm (10x Genomics spec).
# https://kb.10xgenomics.com/s/article/11636252598925-What-are-the-Xenium-image-scale-factors
# Transcript x_location / y_location are in microns.
# To overlay transcripts on the full-res image: pixel = micron / pixel_size.
PIXEL_SIZE_UM = 0.2125


def convert_xenium_to_scs(parquet_path: str,
                          output_bgi_tsv: str,
                          morphology_image_path: str,
                          output_morph2d_tif: str):
    """
    Convert Xenium transcripts to SCS/BGI format with correct pixel-space coordinates.

    Xenium x_location / y_location are in microns.
    The morphology image full resolution is 0.2125 µm/px (Xenium spec).
    Coordinates are converted to pixels: pixel = micron / pixel_size.
    The morphology image is cropped to the pixel ROI covered by the transcripts.
    No pre-binning is applied — binning is handled by SCS itself at runtime.
    """
    transcripts = pd.read_parquet(parquet_path, engine="pyarrow")

    gene_col  = "feature_name"
    x_col     = "x_location"
    y_col     = "y_location"

    table = transcripts[[gene_col, x_col, y_col]].copy()
    table = table.dropna(subset=[gene_col, x_col, y_col])

    # Convert micron coordinates → full-resolution pixel coordinates.
    # Xenium: x_location is along image width (columns), y_location along height (rows).
    # Pixel coordinates are crop-local (zero-based from transcript ROI min).
    r_min = int(round(float(transcripts[y_col].min()) / PIXEL_SIZE_UM))
    c_min = int(round(float(transcripts[x_col].min()) / PIXEL_SIZE_UM))
    table["row"]    = np.rint(table[y_col].astype(float) / PIXEL_SIZE_UM - r_min).astype(int)
    table["column"] = np.rint(table[x_col].astype(float) / PIXEL_SIZE_UM - c_min).astype(int)

    # Each Xenium transcript row is one molecule; aggregate per (gene, pixel).
    table["MIDCounts"] = 1
    table = table.rename(columns={gene_col: "geneID"})
    bgi = (
        table.groupby(["geneID", "row", "column"], as_index=False)["MIDCounts"]
        .sum()
        .rename(columns={"row": "x", "column": "y"})
    )

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

    # Crop to the pixel ROI covered by transcripts (same r_min/c_min as coords above).
    r_max = int(round(float(transcripts[y_col].max()) / PIXEL_SIZE_UM))
    c_max = int(round(float(transcripts[x_col].max()) / PIXEL_SIZE_UM))

    # Clamp to image bounds.
    H, W = image2d.shape
    r_min_clamped = max(0, r_min)
    r_max_clamped = min(H - 1, r_max)
    c_min_clamped = max(0, c_min)
    c_max_clamped = min(W - 1, c_max)

    cropped = image2d[r_min_clamped:r_max_clamped + 1, c_min_clamped:c_max_clamped + 1]

    out_morph2d = Path(output_morph2d_tif)
    out_morph2d.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(out_morph2d, cropped)

    # ── BGI file (SCS/spateo format) ────────────────────────────────────────────
    # spateo read_bgi_agg: x → AnnData dim-0 (height/rows), y → dim-1 (width/cols).
    out_bgi_tsv = Path(output_bgi_tsv)
    out_bgi_tsv.parent.mkdir(parents=True, exist_ok=True)
    bgi.to_csv(out_bgi_tsv, sep="\t", index=False)


if __name__ == "__main__":
    transcripts_parquet: str    = "${transcripts_parquet}"
    morphology_image: str       = "${morphology_image}"
    prefix: str                 = "${prefix}"

    output_bgi_tsv    = f"{prefix}/scs_input_bgi.tsv"
    output_morph2d_tif = f"{prefix}/morph2d.tif"

    convert_xenium_to_scs(
        parquet_path=transcripts_parquet,
        output_bgi_tsv=output_bgi_tsv,
        morphology_image_path=morphology_image,
        output_morph2d_tif=output_morph2d_tif,
    )

    with open("versions.yml", "w", encoding="utf-8") as fobj:
        fobj.write('"${task.process}":\\n')
        fobj.write('xenium2scs: "1.0.0"\\n')
