#!/usr/bin/env python3

import csv
import random
from pathlib import Path


class BaysorPreview():
    """
    Utility class to generate baysor preview dataset
    """
    @staticmethod
    def generate_dataset(
            transcripts: Path,
            sampled_transcripts: Path,
            sample_fraction: float = 0.3,
            random_state: int = 42
        ) -> None:
        """
        Reads a csv file & randomly samples a fraction of rows,
        and writes the result to a .csv file.

        Args:
            transcripts (str): unziped transcripts.csv from xenium bundle
            sampled_transcripts (str): randomly subsampled transcripts.csv file
            sample_fraction (float): Fraction of rows to sample
            random_state (int): Seed for reproducibility
        """

        random.seed(random_state)

        with open(transcripts, mode='rt', newline='') as infile, \
            open(sampled_transcripts, mode='wt', newline='') as outfile:

            reader = csv.reader(infile)
            writer = csv.writer(outfile)

            # get the header libne
            header = next(reader)
            writer.writerow(header)

            # randomize csv rows to write
            for row in reader:
                if random.random() < float(sample_fraction):
                    writer.writerow(row)

        return None

    @staticmethod
    def generate_version_yml() -> None:
        with open("versions.yml", "w") as yml:
            yml.write('"${task.process}":\\n')
            yml.write("Baysor-Preview Create Dataset: 0.7.1'\\n")

        return None


def main() -> None:
    """
    Run create dataset as nf module
    """
    transcripts: str = "${transcripts}"
    sample_fraction: float = "${sample_fraction}"
    sampled_transcripts: str = "sampled_transcripts.csv"

    # generate dataset
    BaysorPreview.generate_dataset (
        transcripts=transcripts,
        sampled_transcripts=sampled_transcripts,
        sample_fraction=sample_fraction
    )

    # generate versions.yml
    BaysorPreview.generate_version_yml()

    return None



if __name__ == "__main__":
    main()
