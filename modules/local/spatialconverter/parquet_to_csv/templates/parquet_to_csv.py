#!/usr/bin/env python

import pandas as pd
from pathlib import Path


def convert_parquet (
        transcripts: Path,
        extension: str = '.csv'
    ) -> None:

    df = pd.read_parquet(transcripts, engine = 'pyarrow')

    if extension == ".gz":
        output = transcripts.replace(".parquet", ".csv.gz")
        df.to_csv(f"{output}", compression='gzip', index=False)
    else:
        output = transcripts.replace(".parquet", ".csv")
        df.to_csv(f"{output}", index=False)

    return None


if __name__ == '__main__':

    transcripts: str = "${transcripts}"
    extension: str = "${extension}"

    # generate transcripts.csv(.gz)
    convert_parquet (
        transcripts=transcripts,
        extension=extension
    )

    #Output versions.yml
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('spatialconverter: "v0.0.1"\\n')
