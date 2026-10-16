"""
Comments to version:
- Uses tissue maps from tissue detector. Therefore, slides should be processed by tissue detector firstly.
- Consider adding color schema if you use the tool for a new entity
"""
from wsi_colors import colors_QC7 as colors
import torch
from PIL import Image
import os
from wsi_slide_info import slide_info
from wsi_process import slide_process_single, mask_to_geojson
from wsi_maps import make_overlay
import numpy as np
import timeit
import cv2
from tqdm import tqdm
import argparse
from pathlib import Path
from wsi_tiff_reader import TiffSlide
Image.MAX_IMAGE_PIXELS = 1000000000

# DEVICE
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
'''
'cuda:0' - NVIDIA GPU card
'mps'    - APPLE Silicon
'''

# Input parameter
parser = argparse.ArgumentParser()
parser.add_argument('--slide_folder', dest='slide_folder', help='path to RGB TIFFs', required=True, type=str)
parser.add_argument('--output_dir', dest='output_dir', help='path to output folder', required=True, type=str)
parser.add_argument('--create_geojson', dest='create_geojson', help='create geojson for QC or not', default="Y", type=str)
parser.add_argument('--start', dest='start', default=0, help='start num of WSIs', type=int)
parser.add_argument('--mpp_model', dest='MPP_MODEL', default=1.5,
                    help='MPP of the training model, should only be 1.0, 1.5, 2.0', type=float)
parser.add_argument('--end', dest='end', default=-1, help='end num of WSIs', type=int)
parser.add_argument('--ol_factor', dest='ol_factor', default=10,
                    help='reduction factor of the overlay compared to dimensions of original WSI', type=int)
parser.add_argument('--mpp', type=float, help='MPP override in microns per pixel; otherwise use OME metadata')
parser.add_argument('--device', choices=('cpu', 'cuda'), default=DEVICE)
parser.add_argument('--model_path', help='path to the artifact model checkpoint')

args = parser.parse_args()

MPP_MODEL = args.MPP_MODEL
start = args.start
end = args.end
SLIDE_DIR = args.slide_folder
OUTPUT_DIR = args.output_dir
OVERLAY_FACTOR = args.ol_factor
create_geojson = args.create_geojson
DEVICE = args.device

# MODEL(S)
# MODEL 1: Artifacts detection

MODEL_QC_DIR = Path(__file__).resolve().parent / 'models/qc'
if MPP_MODEL == 1.5:
    MODEL_QC_NAME = 'GrandQC_MPP15.pth'
elif MPP_MODEL == 1.0:
    MODEL_QC_NAME = 'GrandQC_MPP1.pth'
elif MPP_MODEL == 2.0:
    MODEL_QC_NAME = 'GrandQC_MPP2.pth'
else:
    raise Exception("mpp of the model can only be 1.0, 1.5, 2.0")
M_P_S_MODEL_1 = 512
ENCODER_MODEL_1 = 'timm-efficientnet-b0'
ENCODER_MODEL_1_WEIGHTS = 'imagenet'
BACK_CLASS = 7

if end == -1:
    end = len(os.listdir(SLIDE_DIR))

if create_geojson == "Y":
    geojson_root = os.path.join(OUTPUT_DIR, "geojson_qc")
    os.makedirs(geojson_root, exist_ok=True)
else:
    os.makedirs(OUTPUT_DIR, exist_ok=True)

# OUTPUT (TEXT)
REPORT_FILE_NAME = 'report_slides_' + str(start) + "_" + str(end)   # File name, ".txt" will be added in the end.
REPORT_OUTPUT_DIR = OUTPUT_DIR  # where to save the text report
obj_power = 99

# =============================================================================
# LOAD MODELS
# =============================================================================
model_path = Path(args.model_path) if args.model_path else MODEL_QC_DIR / MODEL_QC_NAME
if not model_path.is_file():
    model_path = (
        Path(__file__).resolve().parent.parent
        / '01_WSI_inference_OPENSLIDE_QC/models/qc'
        / MODEL_QC_NAME
    )
if not model_path.is_file():
    raise FileNotFoundError(f'Artifact detection checkpoint not found: {model_path}')
model = torch.load(model_path, map_location=DEVICE)

# ====================================================================
# PREPARE REPORT FILE, OUTPUT FOLDERS
# =============================================================================

# Prepare report file header
path_result = os.path.join(REPORT_OUTPUT_DIR, REPORT_FILE_NAME + "_stats_per_slide.txt")
output_header = "slide_name" + "\t" + "obj_power" + "\t" + "mpp" + "\t"
output_header = output_header + "patch_n_h_l0" + "\t" + "patch_n_w_l0" + "\t"
output_header = output_header + "patch_overall" + "\t"
output_header = output_header + "width" + "\t" + "height" + "\t"
output_header = output_header + "Time"
output_header = output_header + "\n"
results = open(path_result, "a+")
results.write(output_header)
results.close()

maps_dir = os.path.join(OUTPUT_DIR, 'maps_qc')
overlay_dir = os.path.join(OUTPUT_DIR, 'overlays_qc')
mask_dir = os.path.join(OUTPUT_DIR, 'mask_qc')

os.makedirs(maps_dir, exist_ok=True)
os.makedirs(overlay_dir, exist_ok=True)
os.makedirs(mask_dir, exist_ok=True)

# ====================================================================
# MAIN SCRIPT
# =============================================================================

# Read in flat and OME-TIFF slide names
slide_names = sorted(
    name for name in os.listdir(SLIDE_DIR)
    if name.lower().endswith(('.tif', '.tiff'))
)

# Start analysis loop
for slide_name in slide_names[start:end]:
    # Register start time
    start = timeit.default_timer()
    slide = TiffSlide(os.path.join(SLIDE_DIR, slide_name), args.mpp)
    try:
        p_s, patch_n_w_l0, patch_n_h_l0, mpp, w_l0, h_l0 = slide_info(
            slide, M_P_S_MODEL_1, MPP_MODEL
        )

        tissue_mask_path = os.path.join(
            OUTPUT_DIR, "tis_det_mask", slide_name + '_MASK.png'
        )
        if not os.path.isfile(tissue_mask_path):
            raise FileNotFoundError(
                f"Tissue mask not found: {tissue_mask_path}; run wsi_tis_detect.py first"
            )
        with Image.open(tissue_mask_path) as tis_det_map:
            target_size = (
                max(1, round(w_l0 * mpp / MPP_MODEL)),
                max(1, round(h_l0 * mpp / MPP_MODEL)),
            )
            tis_det_map_mpp = np.asarray(
                tis_det_map.resize(target_size, Image.Resampling.NEAREST)
            )

        map, full_mask = slide_process_single(
            model, tis_det_map_mpp, slide, patch_n_w_l0, patch_n_h_l0, p_s,
            M_P_S_MODEL_1, colors, ENCODER_MODEL_1,
            ENCODER_MODEL_1_WEIGHTS, DEVICE, BACK_CLASS
        )
        stop = timeit.default_timer()
        map_path = os.path.join(maps_dir, slide_name + "_map_ALL.png")
        map.save(map_path)

        mask_path = os.path.join(mask_dir, slide_name + "_mask.png")
        cv2.imwrite(mask_path, full_mask)
        if create_geojson == "Y":
            geojson_path = os.path.join(geojson_root, slide_name + '.geojson')
            factor = MPP_MODEL / mpp
            mask_to_geojson(mask_path, geojson_path, factor)

        overlay = make_overlay(
            slide, map, p_s, patch_n_w_l0, patch_n_h_l0, OVERLAY_FACTOR
        )
        overlay_im = Image.fromarray(overlay)
        overlay_im_name = os.path.join(overlay_dir, slide_name + "_overlay_QC.jpg")
        overlay_im.save(overlay_im_name)

        output_temp = slide_name + "\t"
        output_temp = output_temp + str(obj_power) + "\t" + str(mpp) + "\t"
        output_temp = output_temp + str(patch_n_h_l0) + "\t" + str(patch_n_w_l0) + "\t"
        output_temp = output_temp + str(patch_n_h_l0 * patch_n_w_l0) + "\t"
        output_temp = output_temp + str(patch_n_w_l0 * p_s) + "\t" + str(patch_n_h_l0 * p_s) + "\t"
        output_temp = output_temp + str(round((stop - start) / 60, 1))
        output_temp = output_temp + "\n"

        with open(path_result, "a+") as results:
            results.write(output_temp)
    finally:
        slide.close()
