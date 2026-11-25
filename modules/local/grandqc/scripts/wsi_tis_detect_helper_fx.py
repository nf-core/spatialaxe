import numpy as np

def to_tensor_x(x, **kwargs):
    return x.transpose(2, 0, 1).astype('float32')

def get_preprocessing(image,preprocessing_fn):
    image = np.array(image)
    x = preprocessing_fn(image)
    x = to_tensor_x(x)
    return x

def make_class_map (mask, class_colors):
    r = np.zeros_like(mask).astype(np.uint8)
    g = np.zeros_like(mask).astype(np.uint8)
    b = np.zeros_like(mask).astype(np.uint8)
    for l in range(0, len(class_colors)):
        idx = mask == l
        r[idx] = class_colors [l][0]
        g[idx] = class_colors [l][1]
        b[idx] = class_colors [l][2]

    rgb = np.stack([r, g, b], axis=2)
    return rgb