# MAIN LOOP TO PROCESS WSI
import numpy as np
from PIL import Image
import segmentation_models_pytorch as smp
import torch
from tqdm import tqdm
import cv2
import json


#Helper functions
def to_tensor_x(x, **kwargs):
    return x.transpose(2, 0, 1).astype('float32')

def get_preprocessing(image, preprocessing_fn, model_size):
    if image.size != model_size:
        image = image.resize(model_size)
        print('resized')
    image = np.array(image)
    x = preprocessing_fn(image)
    x = to_tensor_x(x)
    return x

def make_1class_map_thr (mask, class_colors):
    r = np.zeros_like(mask).astype(np.uint8)
    g = np.zeros_like(mask).astype(np.uint8)
    b = np.zeros_like(mask).astype(np.uint8)
    for l in range(1, len(class_colors)+1):
        idx = mask == l
        r[idx] = class_colors [l-1][0]
        g[idx] = class_colors [l-1][1]
        b[idx] = class_colors [l-1][2]

    rgb = np.stack([r, g, b], axis=2)
    return rgb


def slide_process_single(model, tis_det_map_mpp, slide, patch_n_w_l0, patch_n_h_l0, p_s, m_p_s, colors,
                         ENCODER_MODEL_1,ENCODER_WEIGHTS, DEVICE, BACK_CLASS):
    '''
    Tissue detection map is generated under MPP = 4, therefore model patch size of (512,512) corresponds to tis_det_map patch
    size of (128,128).
    '''

    model_size = (m_p_s, m_p_s)
    preprocessing_fn = smp.encoders.get_preprocessing_fn(ENCODER_MODEL_1, ENCODER_WEIGHTS)

    # Start loop
    for he in tqdm(range(patch_n_h_l0), total=patch_n_h_l0):
        h = he * p_s
        # print("Current cycle ", he + 1, " of ", patch_n_h_l0)
        for wi in range(patch_n_w_l0):
            w = wi * p_s
            td_patch = tis_det_map_mpp[he*m_p_s:(he+1)*m_p_s, wi*m_p_s:(wi+1)*m_p_s]
            if td_patch.shape != (512,512):
                # td_patch padding (incase td_patch does not equal (512,512))
                original_shape = td_patch.shape

                # Desired shape
                desired_shape = (512, 512)

                # Calculate padding needed
                padding = [
                    (0, max(0, desired_shape[i] - original_shape[i]))
                    for i in range(2)
                ]

                # Apply padding
                td_patch = np.pad(td_patch, padding, mode='constant')

            if np.count_nonzero(td_patch == 0) > 50:    # here change to check of segmentation map
                # Generate patch
                # Extract the region from the NumPy array
                work_patch_np = slide.read_patch(w, h, p_s, p_s)
                # print(work_patch_np.shape)

                # Create a PIL Image from the NumPy array
                work_patch_pil = Image.new('RGB', (p_s, p_s), 'white')
                work_patch_pil.paste(Image.fromarray(work_patch_np), (0, 0))

                # Resize the image to the model patch size (m_p_s)
                work_patch_resized = work_patch_pil.resize((m_p_s, m_p_s), Image.Resampling.LANCZOS)

                # If needed, convert the image to RGB mode
                work_patch_resized_rgb = work_patch_resized.convert('RGB')

                work_patch = work_patch_resized_rgb

                image_pre = get_preprocessing(work_patch, preprocessing_fn, model_size)
                x_tensor = torch.from_numpy(image_pre).to(DEVICE).unsqueeze(0)
                predictions = model.predict(x_tensor)
                predictions = (predictions.squeeze().cpu().numpy())

                mask_raw = np.argmax(predictions, axis=0).astype('int8')
                mask = np.where(td_patch == 1, BACK_CLASS, mask_raw)


            else:
                mask = np.full((512, 512), BACK_CLASS)

            if wi == 0:
                temp_image = mask

            else:
                temp_image = np.concatenate((temp_image, mask), axis=1)

        if he == 0:
            end_image = temp_image

        else:
            end_image = np.concatenate((end_image, temp_image), axis=0)

    end_image_1class = make_1class_map_thr(end_image, colors)
    end_image_1class = Image.fromarray(end_image_1class)
    end_image_1class = end_image_1class.resize((patch_n_w_l0*200, patch_n_h_l0*200), Image.Resampling.LANCZOS) #what is 50 here?

    return end_image_1class, end_image

def mask_to_geojson(mask_path, output_path, scale_factor=1.0):
    """
    Convert a semantic segmentation mask to GeoJSON with coordinate scaling

    Parameters:
    -----------
    mask_path : str
        Path to the input PNG mask file
    output_path : str
        Path to save the output GeoJSON file
    scale_factor : float, optional, should be: model_mpp / slide_mpp
        Factor to scale coordinates by (default: 1.0)

    Returns:
    --------
    None
    """
    # Define class mapping
    CLASS_MAPPING = {
        1: "Normal Tissue",
        2: "Fold",
        3: "Darkspot & Foreign Object",
        4: "PenMarking",
        5: "Edge & Air Bubble",
        6: "OOF",  # Out of Focus
        7: "Background"
    }

    # Read the mask image
    mask = cv2.imread(mask_path, cv2.IMREAD_UNCHANGED)

    # Dictionary to store features for each class
    features = []

    # Iterate through unique class values (1 to 7)
    for class_value in range(2, 7):
        # Create a binary mask for the current class
        class_mask = (mask == class_value).astype(np.uint8) * 255

        # Find contours
        contours, _ = cv2.findContours(class_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        # Convert contours to GeoJSON features
        for contour in contours:
            # Flatten contour and reshape
            contour_points = contour.reshape(-1, 2)

            # Scale coordinates
            scaled_points = contour_points * scale_factor

            # Skip contours with less than 4 points
            if len(scaled_points) < 4:
                # print(f"Skipping contour with {len(scaled_points)} points for class {class_value}")
                continue

            # Ensure polygon is closed by adding first point at the end if needed
            polygon_points = scaled_points.tolist()
            if not np.array_equal(polygon_points[0], polygon_points[-1]):
                polygon_points.append(polygon_points[0])

            # Create feature for this polygon
            feature = {
                "type": "Feature",
                "properties": {
                    "class_id": int(class_value),
                    "classification": CLASS_MAPPING.get(class_value, "Unknown"),
                    "area": cv2.contourArea(contour) * (scale_factor ** 2)
                },
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [polygon_points]
                }
            }

            features.append(feature)

    # Create GeoJSON structure
    geojson = {
        "type": "FeatureCollection",
        "features": features,
        "metadata": {
            "class_mapping": CLASS_MAPPING,
            "scale_factor": scale_factor
        }
    }

    # Write to file
    with open(output_path, 'w') as f:
        json.dump(geojson, f, indent=2)
