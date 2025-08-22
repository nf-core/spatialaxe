#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from typing import List

# Expected columns in transcripts.parquet
REQUIRED_COLUMNS: List[str] = [
    "transcript_id",
    "cell_id",
    "overlaps_nucleus",
    "feature_name",
    "x_location",
    "y_location",
    "z_location",
    "qv",
    "segger_id"
]

def refine_transcripts(parquet_path: str) -> pd.DataFrame:
    """
    Replace the cell_id column with segger_id
    """
    parquet_file = Path(parquet_path)
    if not parquet_file.exists():
        raise FileNotFoundError(f"File not found: {parquet_path}")

    # Read parquet file
    df = pd.read_parquet(parquet_file, engine="pyarrow")

    # Validate required columns
    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # get 'cell_id' column index
    cell_id_index = df.columns.get_loc("cell_id")

    # Drop 'cell_id' and insert 'segger_id' at the same position
    df = df.drop(columns=["cell_id"])
    segger_series = df.pop("segger_id")
    df.insert(cell_id_index, "cell_id", segger_series)

    return df


def main(input_file: str) -> None:
    transcripts = refine_transcripts(input_file)
    transcripts.to_parquet("transcripts.parquet", engine="pyarrow")


if __name__ == "__main__":

    transcripts = "${transcripts}"

    main(input_file=transcripts)

    #Output versions.yml
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('segger2xr: "v0.0.1"\\n')
