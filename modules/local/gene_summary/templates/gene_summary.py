#!/usr/bin/env python3

# versions
import platform
import yaml

import os

os.environ["MPLCONFIGDIR"] = "./tmp/mpl"
os.environ["NUMBA_CACHE_DIR"] = "./tmp/numba"

import pandas as pd
import scanpy as sc
import numpy as np

from pathlib import Path
from typing import Tuple, List


class Arguments:
    """Parses the arguments, including the ones coming from $task.ext.args.
    Adopted from mygene module (Suzanne Jin)."""

    def __init__(self) -> None:
        self.singlet_str = "singlet"
        self.doublet_str = "doublet"
        self.negative_str = "negative"
        self.parse_input_args()
        self.creat_output_dirs()

    def parse_input_args(self) -> None:
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"

        self.barcodes = "${barcodes}"
        self.vireo = "${vireo}"
        self.demuxlet = "${demuxlet}"
        self.freemuxlet = "${freemuxlet}"
        self.souporcell = "${souporcell}"

        path_vars = {
            "barcodes",
            "vireo",
            "demuxlet",
            "freemuxlet",
            "souporcell",
        }

        def _tranlate_to_python(input_str, value_str):
            if value_str.strip() == "":
                return None
            else:
                path = Path(value_str)
                if not path.exists():
                    raise FileNotFoundError(f"Path does not exist: {path}")
                return path

        for var in path_vars:
            raw_value = getattr(self, var)
            processed_value = _tranlate_to_python(var, raw_value)
            setattr(self, var, processed_value)

    def creat_output_dirs(self) -> None:
        directories = {
            "assignment": "_genetic_summary_assignment.csv",
            "classification": "_genetic_summary_classification.csv",
            "overview_assignment": "_genetic_overview_assignment.csv",
            "overview_classification": "_genetic_overview_classification.csv",
        }

        for output, directory in directories.items():
            setattr(self, output, self.prefix + directory)

    def print_args(self) -> None:
        for attr in vars(self):
            print(f"{attr}: {getattr(self, attr)}")


class ProcessDeconvolutionMethodResult:
    def __init__(self):
        self.deconvolution_methods = ["demuxlet", "freemuxlet", "souporcell", "vireo"]

        self.checkHashNames = True
        self.chechEmptyInput = True

    def vireo(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        results = pd.read_csv(args.vireo, sep="\t")

        assignment = results[["cell", "donor_id"]].rename(
            columns={"cell": "Barcode", "donor_id": "vireo"}
        )
        assignment["vireo"].replace({"unassigned": args.negative_str}, inplace=True)

        classification = assignment.copy()
        classification["vireo"][
            ~classification["vireo"].isin([args.doublet_str, args.negative_str])
        ] = args.singlet_str

        return assignment, classification

    def souporcell(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        results = pd.read_csv(args.souporcell, sep="\t").iloc[:, 0:3]
        results.loc[results["status"] == "doublet", "assignment"] = "doublet"
        results.loc[results["status"] == "unassigned", "assignment"] = "negative"

        assignment = results[["barcode", "assignment"]].rename(
            columns={"barcode": "Barcode", "assignment": "souporcell"}
        )

        classification = assignment.copy()
        classification["souporcell"] = classification["souporcell"].where(
            classification["souporcell"].isin(["doublet", "negative"]), "singlet"
        )

        return assignment, classification

    def demuxlet(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        return self.demuxlet_or_freemuxlet(args, "demuxlet")

    def freemuxlet(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        return self.demuxlet_or_freemuxlet(args, "freemuxlet")

    def demuxlet_or_freemuxlet(
        self, args: Arguments, method: str
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        file_path = getattr(args, method)

        result = pd.read_csv(file_path, sep="\t")

        result[method] = np.where(
            result["BEST.GUESS"].str.split(",").str[0]
            == result["BEST.GUESS"].str.split(",").str[1],
            result["BEST.GUESS"].str.split(",").str[0],
            args.doublet_str,
        )

        result[method] = np.where(
            result["DROPLET.TYPE"] == "AMB", args.negative_str, result[method]
        )

        assignment = result[["BARCODE", method]].rename(columns={"BARCODE": "Barcode"})

        classification = assignment.copy()
        classification[method] = np.where(
            classification[method].isin([args.doublet_str, args.negative_str]),
            classification[method],
            args.singlet_str,
        )

        return assignment, classification


def create_overview_table(dfs: List[pd.DataFrame]):
    """
    Takes the list of assignment/classification DataFrames (assignments/classifications) and prints a summary table:
    method name | total count | match_method1 | match_method2 | ... | count(item1) | count(item2) | ...
    Match to a method counts the number of barcodes that a method has in common with another method.
    An item refers to the donor label in the assignment (HTO-1, HTO-2, ...) or the classification (singlet, doublet, negative).
    """
    rows = []
    all_items = set()
    match_cols = set()

    # extract items and their counts for every deconvolution method
    for df in dfs:
        # add method name and number of barcodes
        method_name = df.columns[1]
        total = len(df)
        row = {"method": method_name, "count_overall": total}

        # add the number of matching barcodes to the other methods
        for df2 in dfs:
            method_name_2 = df2.columns[1]
            match_col_name = f"match_{method_name_2}"
            match_cols.add(match_col_name)
            new_match_col = {
                match_col_name: len(pd.merge(df, df2, on="Barcode", how="inner"))
            }
            row.update(new_match_col)

        # add the counts for each item
        counts = df[method_name].value_counts(dropna=False)
        all_items.update(counts.index)
        row.update(counts.to_dict())

        rows.append(row)

    summary = pd.DataFrame(rows).fillna(0)

    # convert all numeric values to int
    for col in summary.columns:
        if col != "method":
            summary[col] = summary[col].astype(int)

    # order columns
    summary = summary[
        ["method", "count_overall"] + sorted(match_cols) + sorted(list(all_items))
    ]

    return summary


if __name__ == "__main__":
    # ======================== process nextflow input arguments ========================

    args = Arguments()

    # ========================= process results from modules ===========================

    assignments = []
    classifications = []

    # call all functions that process the module outputs
    processing_functions = ProcessDeconvolutionMethodResult()
    for method in processing_functions.deconvolution_methods:
        if getattr(args, method) is not None:
            assignment, classification = getattr(processing_functions, method)(args)
            assignments.append(assignment)
            classifications.append(classification)

    # ================================== save results ==================================

    # ----------------------------------- save csv's -----------------------------------

    # save overview tables
    overview_assignment = create_overview_table(assignments)
    overview_assignment.to_csv(args.overview_assignment, index=False)
    overview_classifications = create_overview_table(classifications)
    overview_classifications.to_csv(args.overview_classification, index=False)

    # save summary of all deconvolution methods
    barcodes_df = pd.read_csv(args.barcodes, header=None, sep="\t", names=["Barcode"])

    # Use barcodes.tsv as index to perform a left join
    assignment_summary = barcodes_df.copy()
    classification_summary = barcodes_df.copy()

    for assignment in assignments:
        assignment_summary = pd.merge(
            assignment_summary, assignment, on="Barcode", how="left"
        )

    for classification in classifications:
        classification_summary = pd.merge(
            classification_summary, classification, on="Barcode", how="left"
        )

    assignment_summary.fillna(args.negative_str).to_csv(args.assignment, index=False)
    classification_summary.fillna(args.negative_str).to_csv(
        args.classification, index=False
    )

    # -------------------------------------- versions ----------------------------------

    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "scanpy": sc.__version__,
            "numpy": np.__version__,
        }
    }

    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)
