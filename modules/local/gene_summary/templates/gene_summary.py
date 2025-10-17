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

def printProccedOutput() -> None:
    # TODO add a function that shows and maybe checks processed results before joining
    print("----- Assignments -----")
    print("")

    for assignment in assignments:
        counts = assignment[assignment.columns[1]].value_counts()
        length = len(assignment)
        print(counts)
        print("length: ", length)
        print("")

    print("----- Classifications -----")
    print("")

    for classification in classifications:
        counts = classification[classification.columns[1]].value_counts()
        length = len(classification)
        print(counts)
        print("length: ", length)
        print("")

if __name__ == "__main__":

    # ======================== process nextflow input arguments ========================

    args = Arguments()
    args.print_args()


    # ========================= process results from modules ===========================

    rna_data = sc.read_10x_mtx(args.rna_matrix)
    hto_data = sc.read_10x_mtx(args.hto_matrix, gex_only=False)

    # call all functions that process the module outptus

    assignments = []
    classifications = []

    processing_functions = ProcessDeconvolutionMethodResult()
    for method in list(processing_functions.deconvolution_methods):
        if getattr(args,method) is not None:
            assignment, classification = getattr(processing_functions, method)(args)
            assignments.append(assignment)
            classifications.append(classification)

    printProccedOutput()

    # ================================== save results ==================================

    # ----------------------------------- save csv's -----------------------------------

    # TODO restructure the if statement if I keep using the hto_data
    # have to to this because demuxem has more barcodes as output that it received as input
    # https://github.com/lilab-bcb/demuxEM/issues/20

    # Read the file as a single-column DataFrame and set the index
    assignment_summary = pd.read_csv(args.barcodes, header=None, names=["Barcode"])
    classification_summary = assignment_summary.copy()

    print(assignment_summary)

    for assignment in assignments:
        assignment_summary = pd.merge(assignment_summary, assignment, on="Barcode", how="left").replace("", args.negative_str).fillna(args.negative_str)

    assignment_summary.to_csv(args.assignment, index=False)

    for classification in classifications:
        classification_summary = pd.merge(classification_summary, classification, on="Barcode", how="left")

    classification_summary.to_csv(args.classification, index=False)

    print(assignment_summary)
    print(classification_summary)
    # -------------------------------- save mudata/anndata -----------------------------

    # if args.generate_mudata or args.generate_anndata:
    #     # join on index (Barcode)
    #     rna_data.obs = rna_data.obs.join(assignment_summary, how="left")
    #     # fill all empty of the used modules with negative values (for expression data)
    #     used_modules = list(assignment_summary.columns)
    #     for col in used_modules:
    #         if pd.api.types.is_categorical_dtype(rna_data.obs[col]):
    #             if args.negative_str not in rna_data.obs[col].cat.categories:
    #                 rna_data.obs[col] = rna_data.obs[col].cat.add_categories([args.negative_str])

    #     rna_data.obs[used_modules] = rna_data.obs[used_modules].fillna(args.negative_str)
    #     rna_data.obs[used_modules] = rna_data.obs[used_modules].astype(str)

    #     if args.generate_mudata:
    #         # join on index (Barcode) and create a mudata object
    #         hto_data.obs = hto_data.obs.join(assignment_summary, how="left")
    #         mudata = MuData({"rna": rna_data, "hto": hto_data})
    #         # TODO mudata update?
    #         mudata.write(args.h5mu)

    #     if args.generate_anndata:
    #         rna_data.write(args.h5ad)

    # -------------------------------------- versions ----------------------------------

    # TODO vervollständigen

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
