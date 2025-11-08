#!/usr/bin/env python3

# versions
import platform
import yaml

import pandas as pd
import scanpy as sc
import numpy as np
import mudata as md
import pegasusio as io

from pathlib import Path
from mudata import MuData
from typing import Tuple

class Arguments:
    # adopted from mygene module (Suzanne Jin)
    """
    Parses the arguments, including the ones coming from $task.ext.args.
    """

    def __init__(self) -> None:

        self.singlet_str = "singlet"
        self.doublet_str = "doublet"
        self.negative_str = "negative"
        self.parse_input_args()
        self.creat_output_dirs()
        #self.testing_inputs()

    def parse_input_args(self) -> None:

        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"

        self.rna_matrix               = "${rna_matrix}"
        self.hto_matrix               = "${hto_matrix}"
        self.barcodes                 = "${barcodes}"
        self.vireo                    = "${vireo}"
        self.demuxlet                 = "${demuxlet}"
        self.freemuxlet               = "${freemuxlet}"
        self.souporcell               = "${souporcell}"

        self.generate_anndata         = "${generate_anndata}"
        self.generate_mudata          = "${generate_mudata}"

        path_vars = {
            "rna_matrix",
            "hto_matrix",
            "barcodes",
            "vireo",
            "demuxlet",
            "freemuxlet",
            "souporcell"
        }

        boolean_vars = {
            "generate_anndata",
            "generate_mudata"
        }

        def _tranlate_to_python(input_str,value_str):
            if value_str.strip() == "":
                return None
            else:
                if input_str in path_vars:
                    return Path(value_str)
                elif input_str in boolean_vars:
                    if value_str == "true":
                        return True
                    else:
                        return False

        vars = path_vars | boolean_vars

        for var in vars:
            raw_value = getattr(self, var)
            processed_value = _tranlate_to_python(var, raw_value)
            setattr(self, var, processed_value)

    def creat_output_dirs(self) -> None:
        directories = {
            'assignment':     '_genetic_summary_assignment.csv',
            'classification': '_genetic_summary_classification.csv',
            'h5mu':           '_genetic_summary.h5mu',
            'h5ad':           '_genetic_summary.h5ad'
        }

        for output, directory in directories.items():
            setattr(self, output, self.prefix + directory)

    # TODO add testing for the inputs
    #def testing_inputs(self) -> None:

    def print_args(self) -> None:
        """
        Print the arguments.
        """
        for attr in vars(self):
            print(f"{attr}: {getattr(self, attr)}")

class ProcessDeconvolutionMethodResult:

    def __init__(self):

        self.deconvolution_methods = {
            "vireo",
            "demuxlet",
            "freemuxlet",
            "souporcell"
        }

        self.checkHashNames = True
        self.chechEmptyInput = True

    def vireo(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:

        results = pd.read_csv(args.vireo, sep="\t")

        assignment = results[['cell','donor_id']].rename(columns={'cell':'Barcode','donor_id': 'vireo'})
        assignment['vireo'].replace({"unassigned": args.negative_str}, inplace=True)

        classification = assignment.copy()
        classification['vireo'][~classification['vireo'].isin([args.doublet_str, args.negative_str])] = args.singlet_str

        print(assignment)
        print(classification)

        return assignment, classification

    def souporcell(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        results = pd.read_csv(args.souporcell, sep="\t").iloc[:, 0:3]
        results.loc[results["status"] == "doublet", "assignment"] = "doublet"
        results.loc[results["status"] == "unassigned", "assignment"] = "negative"

        assignment = results[["barcode", "assignment"]].rename(columns={'barcode':'Barcode', 'assignment': 'souporcell'})

        classification = assignment.copy()
        classification["souporcell"] = classification["souporcell"].where(
            classification["souporcell"].isin(["doublet", "negative"]),
            "singlet"
        )

        print(assignment)
        print(classification)

        return assignment, classification

    def demuxlet(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        return self.demuxlet_or_freemuxlet(args,"demuxlet")

    def freemuxlet(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        return self.demuxlet_or_freemuxlet(args,"freemuxlet")

    def demuxlet_or_freemuxlet(self, args: Arguments, method: str) -> Tuple[pd.DataFrame, pd.DataFrame]:

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

        assignment = result[['BARCODE',method]].rename(columns={'BARCODE':'Barcode'})

        classification = assignment.copy()
        classification[method] = np.where(
            classification[method].isin([args.doublet_str, args.negative_str]),
            classification[method],   # keep original value if in the list
            args.singlet_str                # otherwise set to singlet
        )

        print(assignment)
        print(classification)

        return assignment, classification

#TODO if we keep saving AnnData/MuData in gene/hash_summary add AnnData to container for input type
def saveAnnDataMuData(args: Arguments, assignment_summary: pd.DataFrame, rna_data, hto_data):
    if args.generate_mudata or args.generate_anndata:
        assignment_summary.set_index("Barcode", inplace=True)
        rna_data.obs = rna_data.obs.join(assignment_summary, how="left").fillna(args.negative_str)

    if args.generate_anndata:
        rna_data.write(args.h5ad)

    if args.generate_mudata:
        hto_data.obs = hto_data.obs.join(assignment_summary, how="left")
        mudata = MuData({"rna": rna_data, "hto": hto_data})
        mudata.update()
        mudata.write(args.h5mu)

def print_method_item_counts(dfs):
    """
    Takes the list of assignment/classification DataFrames and prints a summary table:
      method name | total count | count(item1) | count(item2) | ...
    """
    rows = []
    all_items = set()

    for df in dfs:
        # Get second column name
        method_col = df.columns[1]
        # Count occurrences
        counts = df[method_col].value_counts(dropna=False)
        total = len(df)
        all_items.update(counts.index)
        # Build row
        row = {'method': method_col, 'count_overall': total}
        row.update(counts.to_dict())
        rows.append(row)

    # Build dataframe and fill missing item columns
    summary = pd.DataFrame(rows).fillna(0)

    # Convert all numeric values to int
    for col in summary.columns:
        if col != 'method':
            summary[col] = summary[col].astype(int)

    # Order columns
    item_cols = [c for c in summary.columns if c not in ['method', 'count_overall']]
    summary = summary[['method', 'count_overall'] + sorted(item_cols)]

    # Print
    print(summary.to_string(index=False))


if __name__ == "__main__":

    # ======================== process nextflow input arguments ========================

    args = Arguments()

    # only print for debugging
    args.print_args()

    # ========================= process results from modules ===========================

    rna_data = sc.read_10x_mtx(args.rna_matrix)
    hto_data = sc.read_10x_mtx(args.hto_matrix, gex_only=False)

    # call all functions that process the module outputs

    assignments = []
    classifications = []

    processing_functions = ProcessDeconvolutionMethodResult()
    for method in list(processing_functions.deconvolution_methods):
        if getattr(args,method) is not None:
            assignment, classification = getattr(processing_functions, method)(args)
            assignments.append(assignment)
            classifications.append(classification)

    # only print for debugging
    print_method_item_counts(assignments)
    print_method_item_counts(classifications)

    # ================================== save results ==================================

    # ----------------------------------- save csv's -----------------------------------

    # Use rna_data.obs_names() as index to perform a left join
    assignment_summary = pd.DataFrame(rna_data.obs_names, columns=['Barcode'])
    classification_summary = assignment_summary.copy()

    for assignment in assignments:
        assignment_summary = (
            pd.merge(assignment_summary, assignment, on="Barcode", how="left")
              .fillna(args.negative_str)
        )

    for classification in classifications:
        classification_summary = (
            pd.merge(classification_summary, classification, on="Barcode", how="left")
              .fillna(args.negative_str)
        )

    assignment_summary.to_csv(args.assignment, index=False)
    classification_summary.to_csv(args.classification, index=False)

    # -------------------------------- save mudata/anndata -----------------------------

    saveAnnDataMuData(args, assignment_summary,rna_data,hto_data)


    # -------------------------------------- versions ----------------------------------

    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "scanpy": sc.__version__,
            "numpy": np.__version__,
            "mudata": md.__version__,
            "pegasusio": io.__version__,
            "yaml": yaml.__version__,
            }
    }

    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)
