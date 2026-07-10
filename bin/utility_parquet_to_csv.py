#!/usr/bin/env python3
"""
Stream-convert a Parquet file to CSV format using pyarrow batched I/O.

Pure format transformation: every input row maps to one output row.
Memory usage is bounded by --batch-size (default 200_000 rows) rather
than the full row count, so this works on large bundles (e.g. the 10x
Atera WTA panel: 236M rows x 13 cols) within ~200 MB of RAM rather
than the ~30 GB the previous pandas-based eager loader required.
"""

import argparse
import gzip
from pathlib import Path

import pyarrow.csv as pa_csv
import pyarrow.parquet as pq


def stream_parquet_to_csv(
    transcripts: str,
    extension: str,
    prefix: str,
    batch_size: int,
) -> None:
    """Stream a Parquet file to CSV (optionally gzip-compressed)."""
    Path(prefix).mkdir(parents=True, exist_ok=True)
    pf = pq.ParquetFile(transcripts)

    if extension == ".gz":
        out_path = f"{prefix}/" + transcripts.replace(".parquet", ".csv.gz")
        # pyarrow's CSVWriter writes bytes; wrap a gzip stream.
        sink = gzip.open(out_path, "wb")
    else:
        out_path = f"{prefix}/" + transcripts.replace(".parquet", ".csv")
        sink = open(out_path, "wb")

    try:
        with pa_csv.CSVWriter(sink, pf.schema_arrow) as writer:
            for batch in pf.iter_batches(batch_size=batch_size):
                writer.write_batch(batch)
    finally:
        sink.close()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Stream-convert a Parquet file to CSV format."
    )
    parser.add_argument(
        "--transcripts",
        required=True,
        help="Input parquet filename",
    )
    parser.add_argument(
        "--extension",
        default=".csv",
        help="Output extension: '.csv' or '.gz' (default: .csv)",
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Output directory prefix (sample ID)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=200_000,
        help="Rows per batch (default 200000). Memory ~= batch_size * row_size.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    stream_parquet_to_csv(
        transcripts=args.transcripts,
        extension=args.extension,
        prefix=args.prefix,
        batch_size=args.batch_size,
    )
