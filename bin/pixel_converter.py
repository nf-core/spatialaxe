import numpy as np

def to_pixel_coords(xy_world, pixel_width=1.0):
    """
    Approximate translation of world coordinates to pixel coordinates.

    Args:
        xy_world: (N, 2) array of spatial coordinates (in microns).
        pixel_width: Size of each pixel (e.g. 1.0 = 1 micron per pixel).

    Returns:
        xy_pixel: (N, 2) array of integer pixel positions.
    """
    xy_world = np.asarray(xy_world)
    xy_min = xy_world.min(axis=0)
    xy_pixel = ((xy_world - xy_min) / pixel_width).astype(int)
    return xy_pixel
