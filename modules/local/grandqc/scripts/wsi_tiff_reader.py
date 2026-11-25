# additional script, not originally included in grandQC: https://github.com/cpath-ukk/grandqc/tree/main/02_WSI_inference_OME_TIFF_QC
# reads in RGB TIFF/OME-TIFF images in patches and gets the pixel size from OME metadata or an --mpp override

import math
import xml.etree.ElementTree as ET

import cv2
import numpy as np
import tifffile
import zarr


_UNIT_TO_MICRONS = {
    "µm": 1.0,
    "um": 1.0,
    "micrometer": 1.0,
    "micrometre": 1.0,
    "nm": 0.001,
    "nanometer": 0.001,
    "nanometre": 0.001,
    "mm": 1000.0,
    "millimeter": 1000.0,
    "millimetre": 1000.0,
    "cm": 10000.0,
    "centimeter": 10000.0,
    "centimetre": 10000.0,
    "m": 1000000.0,
    "meter": 1000000.0,
    "metre": 1000000.0,
}


def _ome_mpp(ome_metadata):
    if not ome_metadata:
        return None

    root = ET.fromstring(ome_metadata)
    pixels = root.find(".//{*}Pixels")
    if pixels is None:
        return None

    values = []
    for axis in ("X", "Y"):
        value = pixels.get(f"PhysicalSize{axis}")
        if value is None:
            return None
        unit = pixels.get(f"PhysicalSize{axis}Unit", "µm").lower()
        if unit not in _UNIT_TO_MICRONS:
            raise ValueError(f"Unsupported OME physical size unit: {unit}")
        values.append(float(value) * _UNIT_TO_MICRONS[unit])

    if not math.isclose(values[0], values[1], rel_tol=0.01):
        raise ValueError(f"Non-square pixels are not supported (MPP X/Y: {values})")
    return sum(values) / 2


class TiffSlide:
    """Lazy RGB reader for flat and pyramidal TIFF images."""

    def __init__(self, path, mpp=None):
        self._tiff = tifffile.TiffFile(path)
        self._store = None
        try:
            series = self._tiff.series[0]
            level = series.levels[0]
            self.axes = level.axes
            self._store = level.aszarr()
            self._array = zarr.open(self._store, mode="r")
            if isinstance(self._array, zarr.Group):
                multiscales = self._array.attrs.get("multiscales", [])
                datasets = multiscales[0].get("datasets", []) if multiscales else []
                if not datasets or "path" not in datasets[0]:
                    raise ValueError(
                        f"TIFF pyramid has no first-level array in {path}"
                    )
                self._array = self._array[datasets[0]["path"]]
            self._validate_axes()
            self.mpp = mpp if mpp is not None else _ome_mpp(self._tiff.ome_metadata)
            if self.mpp is None:
                raise ValueError(
                    f"No OME physical pixel size found in {path}; supply --mpp in microns per pixel"
                )
            if not math.isfinite(self.mpp) or self.mpp <= 0:
                raise ValueError(f"MPP must be positive and finite, got {self.mpp}")
        except Exception:
            self.close()
            raise

    def _normalized_axes(self):
        axes = list(self.axes)
        shape = list(self._array.shape)
        for index in reversed(range(len(axes))):
            if axes[index] not in "YXSC":
                if shape[index] != 1:
                    raise ValueError(
                        f"Unsupported non-singleton TIFF axis {axes[index]} in {self.axes}"
                    )
                axes.pop(index)
                shape.pop(index)
        return "".join(axes), shape

    def _validate_axes(self):
        axes, shape = self._normalized_axes()
        channel_axis = "S" if "S" in axes else "C" if "C" in axes else None
        if channel_axis is None or len(axes) != 3 or set(axes) != {"Y", "X", channel_axis}:
            raise ValueError(
                f"Expected an RGB TIFF with Y, X and sample/channel axes; got axes {self.axes}"
            )
        if shape[axes.index(channel_axis)] != 3:
            raise ValueError(
                f"Expected exactly 3 RGB samples/channels; got shape {tuple(shape)}"
            )
        self._channel_axis = channel_axis
        self._height = shape[axes.index("Y")]
        self._width = shape[axes.index("X")]

    @property
    def shape(self):
        return (3, self._height, self._width)

    @property
    def width(self):
        return self._width

    @property
    def height(self):
        return self._height

    def read_patch(self, x, y, width, height):
        x0 = max(0, int(x))
        y0 = max(0, int(y))
        x1 = min(self.width, int(x + width))
        y1 = min(self.height, int(y + height))
        if x1 <= x0 or y1 <= y0:
            return np.empty((0, 0, 3), dtype=np.uint8)

        axes, _ = self._normalized_axes()
        slices = tuple(
            slice(y0, y1)
            if axis == "Y"
            else slice(x0, x1)
            if axis == "X"
            else slice(None)
            if axis in "SC"
            else 0
            for axis in self.axes
        )
        image = np.asarray(self._array[slices])
        image = np.moveaxis(
            image,
            (axes.index("Y"), axes.index("X"), axes.index(self._channel_axis)),
            (0, 1, 2),
        )
        if image.dtype != np.uint8:
            if np.issubdtype(image.dtype, np.integer):
                max_value = np.iinfo(image.dtype).max
                image = np.rint(image.astype(np.float32) * (255.0 / max_value))
            else:
                raise ValueError(f"Unsupported TIFF pixel type: {image.dtype}")
            image = np.clip(image, 0, 255).astype(np.uint8)
        return image

    def read_resized(self, width, height, tile_edge=512):
        width = max(1, int(width))
        height = max(1, int(height))
        resized = np.empty((height, width, 3), dtype=np.uint8)
        scale = max(self.width / width, self.height / height)
        output_tile_edge = max(1, min(tile_edge, int(4096 / scale)))

        for y0 in range(0, height, output_tile_edge):
            y1 = min(height, y0 + output_tile_edge)
            source_y0 = round(y0 * self.height / height)
            source_y1 = round(y1 * self.height / height)
            for x0 in range(0, width, output_tile_edge):
                x1 = min(width, x0 + output_tile_edge)
                source_x0 = round(x0 * self.width / width)
                source_x1 = round(x1 * self.width / width)
                patch = self.read_patch(
                    source_x0,
                    source_y0,
                    source_x1 - source_x0,
                    source_y1 - source_y0,
                )
                resized[y0:y1, x0:x1] = cv2.resize(
                    patch, (x1 - x0, y1 - y0), interpolation=cv2.INTER_AREA
                )
        return resized

    def close(self):
        if self._store is not None:
            self._store.close()
            self._store = None
        if self._tiff is not None:
            self._tiff.close()
            self._tiff = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
