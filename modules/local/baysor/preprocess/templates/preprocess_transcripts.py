#!/usr/bin/env python3

import pandas as pd


def filter_transcripts (
    transcripts: str,
    min_qv: float = 20.0,
    min_x: float = 0.0,
    max_x: float = 24000.0,
    min_y: float = 0.0,
    max_y: float = 24000.0
) -> None:
    """
    Filter transcripts based on the specified thresholds

    Args:
    transcripts - path to transcripts parquet
    ----------------------------------- filters --------------------------------------------
    min_qv - minimum Q-Score to pass filtering (default: 20.0)
    min_x  - only keep transcripts whose x-coordinate is greater than specified limit
             if no limit is specified, the default minimum value will be 0.0
    max_x  - only keep transcripts whose x-coordinate is less than specified limit
             if no limit is specified, the default value will retain all
             transcripts since Xenium slide is <24000 microns in x and y (default: 24000.0)
    min_y  - only keep transcripts whose y-coordinate is greater than specified limit
             if no limit is specified, the default minimum value will be 0.0
    max_y  - only keep transcripts whose y-coordinate is less than specified limit
             if no limit is specified, the default value will retain all
             transcripts since Xenium slide is <24000 microns in x and y (default: 24000.0)
    """
    df = pd.read_parquet(transcripts, engine = 'pyarrow')

    # filter transcripts df with thresholds, ignore negative controls
    filtered_df = df[(df["qv"] >= min_qv) &
                                (df["x_location"] >= min_x) &
                                (df["x_location"] <= max_x) &
                                (df["y_location"] >= min_y) &
                                (df["y_location"] <= max_y) &
                                (~df["feature_name"].str.startswith("NegControlProbe_")) &
                                (~df["feature_name"].str.startswith("antisense_")) &
                                (~df["feature_name"].str.startswith("NegControlCodeword_")) &
                                (~df["feature_name"].str.startswith("BLANK_"))]

    # change cell_id of cell-free transcripts from -1 to 0
    neg_cell_row = filtered_df["cell_id"] == -1
    filtered_df.loc[neg_cell_row,"cell_id"] = 0

    # Output filtered transcripts to parquet
    filtered_df.to_parquet(
        '_'.join(["X"+str(min_x)+"-"+str(max_x), "Y"+str(min_y)+"-"+str(max_y), "filtered_transcripts.parquet"]),
        index=False
    )

    return None


def generate_version_yml() -> None:
    with open("versions.yml", "w") as yml:
        yml.write('"${task.process}":\\n')
        yml.write("Baysor-Preprocess Transcripts: 1.0.0'\\n")

    return None


if __name__ == "__main__":

    transcripts: str = "${transcripts}"

    filter_transcripts (
        transcripts=transcripts,
    )

    generate_version_yml()
