#!/usr/bin/env python3

from pathlib import Path
from importlib.metadata import PackageNotFoundError, version

import numpy as np
from tifffile import imread, imwrite

from scportrait.pipeline.project import Project
from scportrait.pipeline.segmentation.workflows import (
	CytosolSegmentationCellpose,
	DAPISegmentationCellpose,
)


def normalize_input_shape(image: np.ndarray) -> np.ndarray:
	image = np.asarray(image)
	image = np.squeeze(image)

	if image.ndim == 2:
		return image[np.newaxis, ...]

	if image.ndim == 3:
		if image.shape[0] <= 4:
			return image
		if image.shape[-1] <= 4:
			return np.moveaxis(image, -1, 0)
		raise ValueError(f"Unable to infer channel axis from image shape {image.shape}")

	if image.ndim == 4:
		image = image[0]
		if image.shape[0] <= 4:
			return image
		if image.shape[-1] <= 4:
			return np.moveaxis(image, -1, 0)
		raise ValueError(f"Unable to infer channel axis from 4D image shape {image.shape}")

	raise ValueError(f"Unsupported morphology image shape {image.shape}")


def extract_label(label_obj) -> np.ndarray:
	if hasattr(label_obj, "scale0"):
		label_obj = label_obj.scale0.image
	elif hasattr(label_obj, "image"):
		label_obj = label_obj.image

	data = getattr(label_obj, "data", label_obj)
	if hasattr(data, "compute"):
		data = data.compute()

	data = np.asarray(data)
	data = np.squeeze(data)
	return data.astype(np.uint32)


def run_segmentation(
	image_path: str,
	nucleus_only: bool,
	nuclei_out: str,
	cells_out: str,
	project_dir: str,
) -> None:
	image = normalize_input_shape(imread(image_path))

	if nucleus_only:
		segmentation_cls = DAPISegmentationCellpose
		config = {
			"DAPISegmentationCellpose": {
				"input_channels": 1,
				"output_masks": 1,
				"cache": ".",
				"chunk_size": 100,
				"nucleus_segmentation": {"model": "nuclei"},
			}
		}
		image = image[:1, :, :]
		channel_names = ["nucleus"]
	else:
		segmentation_cls = CytosolSegmentationCellpose
		if image.shape[0] == 1:
			image = np.repeat(image, 2, axis=0)
		else:
			image = image[:2, :, :]

		config = {
			"CytosolSegmentationCellpose": {
				"input_channels": 2,
				"output_masks": 2,
				"cache": ".",
				"chunk_size": 100,
				"nucleus_segmentation": {"model": "nuclei"},
				"cytosol_segmentation": {"model": "cyto2"},
				"match_masks": True,
				"filter_masks_size": False,
			}
		}
		channel_names = ["nucleus", "cytosol"]

	project = Project(
		project_location=project_dir,
		config_path=config,
		segmentation_f=segmentation_cls,
		overwrite=True,
		debug=False,
	)
	project.load_input_from_array(image, channel_names=channel_names)
	project.segment()

	if "seg_all_nucleus" in project.sdata:
		imwrite(nuclei_out, extract_label(project.sdata["seg_all_nucleus"]))

	if "seg_all_cytosol" in project.sdata:
		imwrite(cells_out, extract_label(project.sdata["seg_all_cytosol"]))


def generate_versions_yml() -> None:
	try:
		scportrait_version = version("scportrait")
	except PackageNotFoundError:
		scportrait_version = "unknown"

	with open("versions.yml", "w", encoding="utf-8") as handle:
		handle.write('"${task.process}":\n')
		handle.write(f"    scportrait: {scportrait_version}\n")


def main() -> None:
	image_path = "${image}"
	prefix = "${prefix}"
	project_dir = "scportrait_project"
	nucleus_only = "${nucleusOnly}".lower() == "true"
	nuclei_out = f"{prefix}/nuclei_labels.tif"
	cells_out = f"{prefix}/cells_labels.tif"

	Path(prefix).mkdir(parents=True, exist_ok=True)
	Path(project_dir).mkdir(parents=True, exist_ok=True)

	run_segmentation(
		image_path=image_path,
		nucleus_only=nucleus_only,
		nuclei_out=nuclei_out,
		cells_out=cells_out,
		project_dir=project_dir,
	)
	generate_versions_yml()


if __name__ == "__main__":
	main()
