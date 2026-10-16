import numpy as np
from PIL import Image
import cv2


# MAKE OVERLAY: HEATMAP ON REDUCED AND CROPPED SLIDE CLON
def make_overlay(slide, wsi_heatmap_im, p_s, patch_n_w_l0, patch_n_h_l0, overlay_factor):
    width = max(1, int(slide.width // overlay_factor))
    height = max(1, int(slide.height // overlay_factor))
    slide_reduced = slide.read_resized(width, height)
    slide_reduced = Image.fromarray(slide_reduced)

    heatmap_width = max(
        1, int(np.ceil(wsi_heatmap_im.width * p_s / (200 * overlay_factor)))
    )
    heatmap_height = max(
        1, int(np.ceil(wsi_heatmap_im.height * p_s / (200 * overlay_factor)))
    )
    heatmap_temp = wsi_heatmap_im.resize(
        (heatmap_width, heatmap_height), Image.Resampling.LANCZOS
    )
    heatmap_crop = heatmap_temp.crop((0, 0, width, height))
    overlay = cv2.addWeighted(
        np.array(slide_reduced), 0.7, np.array(heatmap_crop), 0.3, 0
    )
    return overlay
