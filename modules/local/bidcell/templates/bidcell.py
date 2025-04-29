#!/usr/bin/env python3

# load variables
cpus = '${task.cpus}'
fp_ref = "${fp_ref}"
fp_pos_markers = "${fp_pos_markers}"
fp_neg_markers = "${fp_neg_markers}"
fp_transcripts = "${fp_transcripts}"
fp_dapi = "${fp_dapi}"
config = """${config}"""
outdir = "${prefix}"

# import packages
import os
# ensure cpu limit is set
os.environ["OMP_NUM_THREADS"] = cpus
os.environ["OPENBLAS_NUM_THREADS"] = cpus
os.environ["MKL_NUM_THREADS"] = cpus
os.environ["VECLIB_MAXIMUM_THREADS"] = cpus
os.environ["NUMEXPR_NUM_THREADS"] = cpus

import torch
torch.set_num_threads(int(cpus))

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

######################### pipeline starts #########################
# create output directory
os.makedirs(outdir, exist_ok=True)

# convert parquet files to csv.gz if needed
if fp_transcripts.endswith(".parquet"):
    fp_transcripts = "transcripts.csv.gz"
    with pa.CompressedOutputStream(fp_transcripts, "gzip") as out:
        table_tx = pp.read_table("${fp_transcripts}")
        pc.write_csv(table_tx, out)

# if null, use a dummy ref and set ref markers weights to zero later
if fp_ref == "" and fp_pos_markers == "" and fp_neg_markers == "":
    null_ref = True

    # define toy ref names
    fp_ref = "sc_breast.csv"
    fp_pos_markers = "sc_breast_markers_pos.csv"
    fp_neg_markers = "sc_breast_markers_neg.csv"

    # read from github and write to file
    pd.read_csv("https://raw.githubusercontent.com/SydneyBioX/BIDCell/refs/tags/v1.0.3/data/sc_references/sc_breast.csv", index_col=0).to_csv(fp_ref)
    pd.read_csv("https://raw.githubusercontent.com/SydneyBioX/BIDCell/refs/tags/v1.0.3/data/sc_references/sc_breast_markers_pos.csv", index_col=0).to_csv(fp_pos_markers)
    pd.read_csv("https://raw.githubusercontent.com/SydneyBioX/BIDCell/refs/tags/v1.0.3/data/sc_references/sc_breast_markers_neg.csv", index_col=0).to_csv(fp_neg_markers)
elif fp_ref != "" and fp_pos_markers != "" and fp_pos_markers != "":
    null_ref = False
else:
    raise ValueError("if provided, fp_ref, fp_pos_markers, and fp_neg_marks must be provided all together.")

# write config to file
with open(f"{outdir}/bidcell_config.yaml", "w") as f:
    f.write(
        f"""# Define the config file
        cpus: {cpus} # number of CPUs for multiprocessing

        files: # NOTE: please ensure these point to the right locations
            data_dir: {outdir} # data directory for processed/output data
            fp_dapi: {fp_dapi} # path of DAPI image or path of output stitched DAPI if using stitch_nuclei
            fp_transcripts: {fp_transcripts} # path of transcripts file
            fp_ref: {fp_ref} # file path of reference data
            fp_pos_markers: {fp_pos_markers} # file path of positive markers
            fp_neg_markers: {fp_neg_markers} # file path of negative markers

        {config}
        """
    )

model = bidcell.BIDCellModel(f"{outdir}/bidcell_config.yaml")

# if null ref, we set marker weights to zero
if null_ref:
    model.config.training_params.pos_weight = 0
    model.config.training_params.neg_weight = 0

if model.config.model_params.elongated == ["placeholder"]:
    no_elongated = True
else:
    no_elongated = False

# process transcripts
# find discovered genes
observed_genes = pd.read_csv(fp_transcripts, usecols=['feature_name'])['feature_name'].unique().tolist()
# if we have pattern to ignore, we remove these genes
if model.config.transcripts.transcripts_to_filter is not None:
    observed_genes = [g for g in observed_genes if not g.startswith(tuple(model.config.transcripts.transcripts_to_filter)) ]

# load reference to compare
df_ref_orig = pd.read_csv(model.config.files.fp_ref, index_col=0)
ct_columns = df_ref_orig.columns[-3:].tolist()

# get genes not in ref
extra_genes = [g for g in observed_genes if g not in df_ref_orig.columns]

# - If we have genes in data but not in ref, we add these to ref with 0 counts,
# - If we do not have elongated cells types, we add a dummy cell type "placeholder" to the ref with 0 counts.
# we write a new ref to replace the original ref if any of these conditions are met
if len(extra_genes) > 0 or no_elongated:
    # define new file paths
    fp_ref = "expanded_ref.csv.gz"
    fp_pos_markers = "expanded_pos_markes.csv.gz"
    fp_neg_markers = "expanded_neg_markes.csv.gz"

    # read in positive and negative markers
    df_ref = pd.read_csv(model.config.files.fp_ref, index_col=0)
    df_pos_markers = pd.read_csv(model.config.files.fp_pos_markers, index_col=0)
    df_neg_markers = pd.read_csv(model.config.files.fp_neg_markers, index_col=0)

    # add extra genes to ref
    if len(extra_genes) > 0:
        # keep only expressed genes and save to files
        df_ref[extra_genes] = 0
        df_ref = df_ref[observed_genes + ct_columns]

        df_pos_markers[extra_genes] = 0
        df_pos_markers = df_pos_markers[observed_genes]

        df_neg_markers[extra_genes] = 0
        df_neg_markers = df_neg_markers[observed_genes]
        print("Found genes in data but not in ref. Added these to ref with zero counts.")

    # If the no_elongated flag is set, we need to add a row called "placeholder" to the three files with all zeros
    if no_elongated:
        # add a row of zeros to the ref
        # df_ref.loc["placeholder"] = 0
        # df_ref.loc["placeholder", ct_columns] = [314159,"placeholder","atlas"]

        # add a row of zeros to the pos markers
        df_pos_markers.loc["placeholder"] = 0

        # add a row of zeros to the neg markers
        df_neg_markers.loc["placeholder"] = 0
        print("No elongated cell types found. Added a dummy cell type 'placeholder' to the ref with zero counts.")

    # write to new files
    df_ref.to_csv(fp_ref, index=True, mode='wb')
    df_pos_markers.to_csv(fp_pos_markers, index=True, mode='wb')
    df_neg_markers.to_csv(fp_neg_markers, index=True, mode='wb')

    # update config paths
    model.config.files.fp_ref = fp_ref
    model.config.files.fp_pos_markers = fp_pos_markers
    model.config.files.fp_neg_markers = fp_neg_markers

# for better control, we run each step separately
# https://github.com/SydneyBioX/BIDCell/blob/e565988cd2e78e622c68bd0a5649a1ec8b9b281f/bidcell/BIDCellModel.py#L35
# preprocess
print("############### Preprocessing ###############")
model.preprocess()
torch.cuda.empty_cache()
print()
print("############### Training ###############")
print()
model.train()
torch.cuda.empty_cache()
print()
print("############### Predict ###############")
torch.cuda.empty_cache()
print()
model.predict()
torch.cuda.empty_cache()
print()
print("############### Done ###############")

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
