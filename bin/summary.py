#!/usr/bin/env python
import os
import argparse
import pandas as pd

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


def merge_dataframes(dataframes):
    merged_df = pd.DataFrame()
    for df in dataframes:
        if merged_df.empty:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on="Barcode", how="outer")
    return merged_df


if __name__ == "__main__":
    if not os.path.exists("summary"):
        os.makedirs("summary")

    # Assignments
    assignment_gene = [
        os.path.join(args.gene_demulti, gene_file)
        for gene_file in os.listdir(args.gene_demulti)
        if gene_file.endswith("_assignment_all.csv") and not gene_file.startswith(".")
    ][0]
    assignment_gene = pd.read_csv(assignment_gene, dtype=str)
    assignment_hash = [
        os.path.join(args.hash_demulti, hash_file)
        for hash_file in os.listdir(args.hash_demulti)
        if hash_file.endswith("_assignment_all.csv") and not hash_file.startswith(".")
    ][0]
    assignment_hash = pd.read_csv(assignment_hash, dtype=str)
    assignment_all = merge_dataframes([assignment_gene, assignment_hash])

    assignment_all = assignment_all.replace({"DBL": "doublet", "AMB": "negative"})
    assignment_all.to_csv(
        "summary/assignment_all_genetic_and_hash.csv", index=False, sep="\t"
    )

    # Classifications
    classification_gene = [
        os.path.join(args.gene_demulti, gene_file)
        for gene_file in os.listdir(args.gene_demulti)
        if gene_file.endswith("_classification_all.csv")
        and not gene_file.startswith(".")
    ][0]
    classification_gene = pd.read_csv(classification_gene, dtype=str)
    classification_hash = [
        os.path.join(args.hash_demulti, hash_file)
        for hash_file in os.listdir(args.hash_demulti)
        if hash_file.endswith("_classification_all.csv")
        and not hash_file.startswith(".")
    ][0]
    classification_hash = pd.read_csv(classification_hash, dtype=str)
    classification_all = merge_dataframes([classification_gene, classification_hash])
    classification_all = classification_all.replace(
        {"SNG": "singlet", "DBL": "doublet", "AMB": "negative"}
    )
    classification_all.to_csv(
        "summary/classification_all_genetic_and_hash.csv", index=False, sep="\t"
    )
