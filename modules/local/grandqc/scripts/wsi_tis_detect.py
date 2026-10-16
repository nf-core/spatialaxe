import os
import numpy as np
import torch
from PIL import Image
import segmentation_models_pytorch as smp
from wsi_tis_detect_helper_fx import to_tensor_x, get_preprocessing, make_class_map
import argparse
from pathlib import Path
from wsi_tiff_reader import TiffSlide
import cv2

# DEVICE
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
'''
'cuda' - NVIDIA GPU card
'mps'    - APPLE Silicon
'''

# MODEL(S)
# MODEL TISSUE DETECTION:
MODEL_TD_DIR = './models/td/'
MODEL_TD_NAME = 'Tissue_Detection_MPP10.pth'
MPP_MODEL_TD = 10
M_P_S_MODEL_TD = 512
ENCODER_MODEL_TD = 'timm-efficientnet-b0'
ENCODER_MODEL_TD_WEIGHTS = 'imagenet'


# OVERLAY PARAMETERS (TRANSPARENCY)
OVER_IMAGE = 0.7    # % original image
OVER_MASK = 0.3     # % segmentation mask

# COLORS for MASK
colors = [[50, 50, 250],    # BLUE: TISSUE
          [128, 128, 128]]  # GRAY: BACKGROUND

parser = argparse.ArgumentParser()
parser.add_argument('--slide_folder', dest='slide_folder', help='path to RGB TIFFs', required=True, type=str)
parser.add_argument('--output_dir', dest='output_dir', help='path to output folder', required=True, type=str)
parser.add_argument('--mpp', dest='mpp', help='MPP override in microns per pixel; otherwise use OME metadata', type=float)
parser.add_argument('--device', choices=('cpu', 'cuda'), default=DEVICE)
parser.add_argument('--model_path', help='path to the tissue model checkpoint')
args = parser.parse_args()

SLIDE_DIR = args.slide_folder
OUTPUT_DIR = args.output_dir
DEVICE = args.device
MODEL_TD_PATH = args.model_path or str(
    Path(__file__).resolve().parent.parent
    / '01_WSI_inference_OPENSLIDE_QC/models/td/Tissue_Detection_MPP10.pth'
)

# Create output dirs
tis_det_dir_mask = os.path.join(OUTPUT_DIR, 'tis_det_mask/')
tis_det_dir_over = os.path.join(OUTPUT_DIR, 'tis_det_overlay/')
tis_det_dir_thumb = os.path.join(OUTPUT_DIR, 'tis_det_thumbnail/')
tis_det_dir_mask_col = os.path.join(OUTPUT_DIR, 'tis_det_mask_col/')

try:
    os.makedirs(tis_det_dir_mask)
    os.makedirs(tis_det_dir_over)
    os.makedirs(tis_det_dir_thumb)
    os.makedirs(tis_det_dir_mask_col)
except:
    print('The folders are already there ..')

# Get slide names
slide_names = sorted(
    name for name in os.listdir(SLIDE_DIR)
    if name.lower().endswith(('.tif', '.tiff'))
)

preprocessing_fn = smp.encoders.get_preprocessing_fn(ENCODER_MODEL_TD, ENCODER_MODEL_TD_WEIGHTS)
model = smp.UnetPlusPlus(
    encoder_name=ENCODER_MODEL_TD,
    encoder_weights=ENCODER_MODEL_TD_WEIGHTS,
    classes=2,
    activation=None,
)
if not os.path.isfile(MODEL_TD_PATH):
    raise FileNotFoundError(f'Tissue detection checkpoint not found: {MODEL_TD_PATH}')
model.load_state_dict(torch.load(MODEL_TD_PATH, map_location='cpu'))
model.to(DEVICE)
model.eval()

# Start analysis loop
for slide_name in slide_names:
    slide = None
    try:
        print("")
        print("Working with: ", slide_name)
        path_slide = os.path.join(SLIDE_DIR, slide_name)
        slide = TiffSlide(path_slide, args.mpp)
        width = max(1, round(slide.width * slide.mpp / MPP_MODEL_TD))
        height = max(1, round(slide.height * slide.mpp / MPP_MODEL_TD))
        image_or = slide.read_resized(width, height)
        Image.fromarray(image_or).save(os.path.join(tis_det_dir_thumb, slide_name + ".jpg"), quality=80)

        '''
        As tissue detector was trained on jpeg compressed images - we have to reproduce this step.
        Otherwise it functions suboptimal.
        '''

        image = image_or
        encode_param = [int(cv2.IMWRITE_JPEG_QUALITY), 80]
        result, image = cv2.imencode('.jpg', image, encode_param)
        image = cv2.imdecode(image, 1)
        image = Image.fromarray(image)

        width, height = image.size

        mask_canvas = np.ones((height, width), dtype=np.uint8)
        class_canvas = np.empty((height, width, 3), dtype=np.uint8)
        class_canvas[:] = colors[1]
        for y in range(0, height, M_P_S_MODEL_TD):
            for x in range(0, width, M_P_S_MODEL_TD):
                image_work = Image.new('RGB', (M_P_S_MODEL_TD, M_P_S_MODEL_TD), 'black')
                image_work.paste(
                    image.crop((x, y, x + M_P_S_MODEL_TD, y + M_P_S_MODEL_TD)),
                    (0, 0),
                )
                image_pre = get_preprocessing(image_work, preprocessing_fn)
                x_tensor = torch.from_numpy(image_pre).to(DEVICE).unsqueeze(0)
                predictions = model.predict(x_tensor)
                predictions = (predictions.squeeze().cpu().numpy())

                mask = np.argmax(predictions, axis=0).astype('int8')
                class_mask = make_class_map(mask, colors)
                patch_height = min(M_P_S_MODEL_TD, height - y)
                patch_width = min(M_P_S_MODEL_TD, width - x)
                mask_canvas[y:y + patch_height, x:x + patch_width] = mask[:patch_height, :patch_width]
                class_canvas[y:y + patch_height, x:x + patch_width] = class_mask[:patch_height, :patch_width]

        Image.fromarray(mask_canvas).save(os.path.join(tis_det_dir_mask, slide_name + '_MASK.png'))
        Image.fromarray(class_canvas).save(os.path.join(tis_det_dir_mask_col, slide_name + '_MASK_COL.png'))
        overlay = cv2.addWeighted(np.array(image), OVER_IMAGE, class_canvas, OVER_MASK, 0)
        overlay = Image.fromarray(overlay)
        overlay.save(os.path.join(tis_det_dir_over, slide_name + '_OVERLAY.jpg'))
    finally:
        if slide is not None:
            slide.close()
