# Adapted GrandQC flat-TIFF scripts

These are copies of the locally adapted GrandQC OME-TIFF workflow scripts. They read single-level RGB TIFF/OME-TIFF files in patches and use OME physical pixel-size metadata or a required `--mpp` override.

Run from this folder so the helper imports resolve:

```bash
cd /workspaces/spatialaxe/modules/local/grandqc/scripts

SLIDES=/path/to/rgb-tiff-folder
OUTPUT=/path/to/grandqc-output
MPP=<measured-microns-per-pixel>

python wsi_tis_detect.py \
  --slide_folder "$SLIDES" \
  --output_dir "$OUTPUT" \
  --mpp "$MPP" \
  --device cpu \
  --model_path /workspaces/spatialaxe/data/grandqc/01_WSI_inference_OPENSLIDE_QC/models/td/Tissue_Detection_MPP10.pth

python main.py \
  --slide_folder "$SLIDES" \
  --output_dir "$OUTPUT" \
  --mpp "$MPP" \
  --mpp_model 1.5 \
  --create_geojson N \
  --device cpu \
  --model_path /workspaces/spatialaxe/data/grandqc/01_WSI_inference_OPENSLIDE_QC/models/qc/GrandQC_MPP15.pth
```

Replace `MPP` with the image's actual physical pixel size. Do not use a guessed value; the web-map H&E mosaic in `data/xenium_test_images` has no calibrated MPP and is only suitable for testing TIFF reading.

Python dependencies are listed in `/workspaces/spatialaxe/data/grandqc/requirements.txt`. These scripts and the GrandQC model checkpoints are subject to GrandQC's non-commercial license; see `/workspaces/spatialaxe/data/grandqc/LICENSE`.
