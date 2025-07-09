#!/usr/bin/env python
import argparse
import pandas as pd
from pathlib import Path

parser = argparse.ArgumentParser(description="Parameters for summarizing results")
parser.add_argument(
    "--gene_demulti",
    help="Folder containing output files of genetic demultiplexing pipeline",
    default=None,
)
parser.add_argument(
    "--hash_demulti",
    help="Folder containing output files of hashing demultiplexing pipeline",
    default=None,
)
args = parser.parse_args()


def merge_dataframes(dataframes: list[pd.DataFrame]) -> pd.DataFrame:
    merged_df = pd.DataFrame()
    for df in dataframes:
        if merged_df.empty:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on="Barcode", how="outer")
    return merged_df


def find_first_file(directory: Path, suffix: str) -> Path:
    return [
        file
        for file in directory.iterdir()
        if file.name.endswith(suffix) and not file.name.startswith(".")
    ][0]


def process_file_pair(
    gene_dir: Path,
    hash_dir: Path,
    suffix: str,
    output_path: Path,
    replacements: dict[str, str],
) -> None:
    gene_file = find_first_file(gene_dir, suffix)
    hash_file = find_first_file(hash_dir, suffix)

    gene_df = pd.read_csv(gene_file, dtype=str)
    hash_df = pd.read_csv(hash_file, dtype=str)

    merged_df = merge_dataframes([gene_df, hash_df])
    merged_df = merged_df.replace(replacements)
    merged_df.to_csv(output_path, index=False, sep="\t")


if __name__ == "__main__":
    summary_dir = Path("summary")
    summary_dir.mkdir(exist_ok=True)

    gene_dir = Path(args.gene_demulti)
    hash_dir = Path(args.hash_demulti)

    process_file_pair(
        gene_dir,
        hash_dir,
        "_assignment_all.csv",
        summary_dir / "assignment_all_genetic_and_hash.csv",
        {"DBL": "doublet", "AMB": "negative"},
    )

    process_file_pair(
        gene_dir,
        hash_dir,
        "_classification_all.csv",
        summary_dir / "classification_all_genetic_and_hash.csv",
        {"SNG": "singlet", "DBL": "doublet", "AMB": "negative"},
    )
