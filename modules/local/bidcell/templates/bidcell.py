#!/usr/bin/env python3

# import packages
import os
import platform
import pyarrow as pa
import pyarrow.parquet as pp
import pyarrow.csv as pc
import bidcell
import pandas as pd
from pip._vendor import pkg_resources

def get_version(package):
    package = package.lower()
    return next((p.version for p in pkg_resources.working_set if p.project_name.lower() == package), "No match")

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.

    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.

    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "    " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

# create output directory
outdir = "${prefix}"
if not os.path.exists(outdir):
    os.makedirs(outdir)

# convert parquet files to csv.gz
fp_transcripts = "${transcripts}"
if "${transcripts}".endswith(".parquet"):
    fp_transcripts = "transcripts.csv.gz"
    with pa.CompressedOutputStream(fp_transcripts, "gzip") as out:
        table_tx = pp.read_table("${transcripts}")
        pc.write_csv(table_tx, out)

# write config to file
with open(f"{outdir}/bidcell_config.yaml", "w") as f:
    f.write(
        f"""# Define the config file
        cpus: ${task.cpus} # number of CPUs for multiprocessing

        files: # NOTE: please ensure these point to the right locations
            data_dir: ${prefix} # data directory for processed/output data
            fp_dapi: ${dapi} # path of DAPI image or path of output stitched DAPI if using stitch_nuclei
            fp_transcripts: {fp_transcripts} # path of transcripts file
            fp_ref: ${ref} # file path of reference data
            fp_pos_markers: ${pos_markers} # file path of positive markers
            fp_neg_markers: ${neg_markers} # file path of negative markers

        ${config}
        """
    )

model = bidcell.BIDCellModel(f"{outdir}/bidcell_config.yaml")

# if reference data is not provided, we fool the model with a dummy one
if "${params.containsKey('pos_markers')}" == "false": 
    model.config.training_params.pos_weight = 0

if "${params.containsKey('neg_markers')}" == "false":
    model.config.training_params.neg_weight = 0

# if elongated is not provided (with the default value placeholder present), we fill it.
if model.config.model_params.elongated == ["placeholder"]:
    model.config.model_params.elongated = pd.read_csv("${ref}").cell_type.unique().tolist()

# run the pipeline
model.run_pipeline()

# write versions to file
versions = {
    "${task.process}": {
        "python": platform.python_version(), 
        "bidcell": get_version("bidcell"), 
        "pyarrow": pa.__version__
    }
}

with open("versions.yml", "w") as f:
    f.write(format_yaml_like(versions))