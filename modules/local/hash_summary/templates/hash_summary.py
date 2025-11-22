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
import mudata as md
import pegasusio as io

from pathlib import Path
from mudata import MuData
from typing import Tuple


class Arguments:
    """Parses the arguments, including the ones coming from $task.ext.args.
    Adopted from mygene module (Suzanne Jin)."""

    def __init__(self) -> None:
        self.singlet_str = "singlet"
        self.doublet_str = "doublet"
        self.negative_str = "negative"
        self.parse_input_args()
        self.creat_output_dirs()
        self.testing_inputs()

    def parse_input_args(self) -> None:
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"

        self.rna_matrix = "${rna_matrix}"
        self.hto_matrix = "${hto_matrix}"
        self.htodemux_assignments = "${htodemux_assignments}"
        self.htodemux_classification = "${htodemux_classification}"
        self.multiseq = "${multiseq}"
        self.bff = "${bff}"
        self.demuxem = "${demuxem}"
        self.gmmdemux_results = "${gmmdemux_results}"
        self.gmmdemux_config = "${gmmdemux_config}"
        self.hasheddrops_results = "${hasheddrops_results}"
        self.hasheddrops_id_to_hash = "${hasheddrops_id_to_hash}"
        self.hashsolo = "${hashsolo}"

        self.generate_anndata = "${generate_anndata}"
        self.generate_mudata = "${generate_mudata}"
        self.bff_methods = "${bff_methods}"
        self.hash_list = "${hash_list}"

        path_vars = {
            "rna_matrix",
            "hto_matrix",
            "htodemux_assignments",
            "htodemux_classification",
            "multiseq",
            "bff",
            "demuxem",
            "gmmdemux_results",
            "gmmdemux_config",
            "hasheddrops_results",
            "hasheddrops_id_to_hash",
            "hashsolo",
        }

        boolean_vars = {"generate_anndata", "generate_mudata"}

        other_vars = {"bff_methods", "hash_list"}

        def _tranlate_to_python(input_str, value_str):
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
                elif input_str == "bff_methods":
                    if value_str == "RAW":
                        return ["bff_raw"]
                    elif value_str == "CLUSTER":
                        return ["bff_cluster"]
                    elif value_str == "COMBINED":
                        return ["bff_raw", "bff_cluster", "bff_consensuscall"]
                    else:
                        raise ValueError(
                            f"Methods ({value_str}) for bff not specified correctly. Choose RAW, CLUSTER or COMBINED as input."
                        )
                elif input_str == "hash_list":
                    return set(
                        hash.strip() for hash in "${hash_list}".strip("[]").split(",")
                    )

        vars = path_vars | boolean_vars | other_vars

        for var in vars:
            raw_value = getattr(self, var)
            processed_value = _tranlate_to_python(var, raw_value)
            setattr(self, var, processed_value)

    def creat_output_dirs(self) -> None:
        directories = {
            "assignment": "_hashing_summary_assignment.csv",
            "classification": "_hashing_summary_classification.csv",
            "h5mu": "_hashing_summary.h5mu",
            "h5ad": "_hashing_summary.h5ad",
        }

        for output, directory in directories.items():
            setattr(self, output, self.prefix + directory)

    def testing_inputs(self) -> None:
        if [self.htodemux_assignments, self.htodemux_classification].count(None) == 1:
            raise ValueError(
                "The assignment or classification file of htodemux is empty."
            )

        if [self.gmmdemux_results, self.gmmdemux_config].count(None) == 1:
            raise ValueError("The results or config file of gmmdemux is empty.")

    def print_args(self) -> None:
        for attr in vars(self):
            print(f"{attr}: {getattr(self, attr)}")


class ProcessModuleOutput:
    def __init__(self):
        # necessary to verify which functions should be called
        # because gmmdemux, hasheddrops and htodemux need two input files
        self.function_name_to_args_name = {
            "demuxem": "demuxem",
            "hashsolo": "hashsolo",
            "hasheddrops": "hasheddrops_results",
            "multiseq": "multiseq",
            "htodemux": "htodemux_assignments",
            "gmmdemux": "gmmdemux_results",
            "bff": "bff",
        }

        self.checkHashNames = True
        self.chechEmptyInput = True

    def demuxem(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        data = io.read_input(str(args.demuxem))
        classification = data.obs["demux_type"].to_frame()
        classification.reset_index(inplace=True)
        classification.columns = ["Barcode", "demuxem"]
        classification["demuxem"] = classification["demuxem"].cat.rename_categories(
            {"unknown": args.negative_str}
        )

        # TODO demuxem has more output barcodes than input barcodes metioned here: https://github.com/lilab-bcb/demuxEM/issues/20
        assignment = data.obs["assignment"].to_frame()
        assignment.reset_index(inplace=True)
        assignment.columns = ["Barcode", "demuxem"]

        return assignment, classification

    def hashsolo(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        results = pd.read_csv(args.hashsolo, index_col=0)
        assignment = results[["Classification"]]

        assignment.columns = ["hashsolo"]
        assignment = assignment.replace(
            {"Doublet": args.doublet_str, "Negative": args.negative_str}
        )

        classification = results[["most_likely_hypothesis"]].copy()
        classification["most_likely_hypothesis"] = classification[
            "most_likely_hypothesis"
        ].replace(
            {0.0: args.negative_str, 1.0: args.singlet_str, 2.0: args.doublet_str}
        )

        classification = classification.rename(
            columns={"most_likely_hypothesis": "hashsolo"}
        )

        return assignment.reset_index(), classification.reset_index()

    def hasheddrops(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        idx_to_htoname_df = pd.read_csv(args.hasheddrops_id_to_hash)
        idx_to_htoname_df.loc[len(idx_to_htoname_df)] = [np.nan, args.negative_str]
        idx_to_htoname_map = idx_to_htoname_df.set_index("Index")["HTO"].to_dict()

        obs_res = pd.read_csv(args.hasheddrops_results)

        obs_res["Classification"] = np.where(
            obs_res["Confident"] & obs_res["Confident"].notna(),
            args.singlet_str,
            np.where(
                obs_res["Doublet"] & obs_res["Doublet"].notna(),
                args.doublet_str,
                args.negative_str,
            ),
        )

        obs_res["Assignment"] = np.where(
            obs_res["Classification"].isin([args.doublet_str, args.negative_str]),
            obs_res["Classification"],
            obs_res["Best"].map(idx_to_htoname_map),
        )

        obs_res.rename(columns={obs_res.columns[0]: "Barcode"}, inplace=True)

        classification = obs_res[["Barcode", "Classification"]].rename(
            columns={"Classification": "hasheddrops"}
        )
        assignment = obs_res[["Barcode", "Assignment"]].rename(
            columns={"Assignment": "hasheddrops"}
        )

        return assignment, classification

    def multiseq(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        assignment = pd.read_csv(args.multiseq)
        assignment.columns = ["Barcode", "multiseq"]
        assignment.replace(
            {"Doublet": args.doublet_str, "Negative": args.negative_str}, inplace=True
        )

        classification = assignment.copy()
        classification.loc[
            (classification["multiseq"] != args.doublet_str)
            & (classification["multiseq"] != args.negative_str),
            "multiseq",
        ] = args.singlet_str

        return assignment, classification

    def htodemux(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        assignment = pd.read_csv(args.htodemux_assignments)
        assignment.columns = ["Barcode", "htodemux"]
        assignment.replace("Doublet", args.doublet_str, inplace=True)
        assignment.replace(
            {"Doublet": args.doublet_str, "Negative": args.negative_str}, inplace=True
        )

        classification = pd.read_csv(args.htodemux_classification)
        classification.columns = ["Barcode", "htodemux"]
        classification.replace(
            {
                "Singlet": args.singlet_str,
                "Doublet": args.doublet_str,
                "Negative": args.negative_str,
            },
            inplace=True,
        )

        return assignment, classification

    def gmmdemux(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        number_of_hashes = len(args.hash_list)

        df_config = pd.read_csv(
            args.gmmdemux_config, header=None, skipinitialspace=True
        )
        df_config.columns = ["Cluster_id", "Description"]

        def _classify_hash(cluster_id: int, number_hashes: int) -> str:
            if cluster_id == 0:
                return args.negative_str
            elif 1 <= cluster_id <= number_hashes:
                return args.singlet_str
            else:
                return args.doublet_str

        df_config["Classification"] = df_config["Cluster_id"].apply(
            lambda cluster_id: _classify_hash(cluster_id, number_of_hashes)
        )

        df_config["Assignment"] = df_config["Description"].where(
            df_config["Classification"] == args.singlet_str,
            other=df_config["Classification"],
        )

        # results with Cluster_id's
        df_results = pd.read_csv(args.gmmdemux_results)
        df_results.columns = ["Barcode", "Cluster_id", "Confidence"]

        df_results = df_results.merge(df_config, on="Cluster_id", how="left")

        assignment = df_results[["Barcode", "Assignment"]]
        assignment.columns = ["Barcode", "gmmdemux"]

        classification = df_results[["Barcode", "Classification"]]
        classification.columns = ["Barcode", "gmmdemux"]

        return assignment, classification

    def bff(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:
        df_result = pd.read_csv(args.bff)
        df_result.rename(columns={"cellbarcode": "Barcode"}, inplace=True)

        if len(args.bff_methods) == 3:
            cols = ["bff_raw", "bff_cluster", "consensuscall", "consensuscall.global"]
        else:
            cols = args.bff_methods

        df_result[cols] = df_result[cols].replace(
            {
                "Singlet": args.singlet_str,
                "Doublet": args.doublet_str,
                "Negative": args.negative_str,
                "Discordant": args.negative_str,
                "Not Called": args.negative_str,
            }
        )

        if len(args.bff_methods) == 3:
            # use the classification of consensuscall.global
            assignment = df_result[
                ["Barcode", "bff_raw", "bff_cluster", "consensuscall"]
            ].rename(columns={"consensuscall": "bff_consensuscall"})

            classification = df_result[
                ["Barcode", "bff_raw", "bff_cluster", "consensuscall.global"]
            ].rename(columns={"consensuscall.global": "bff_consensuscall"})
        else:
            assignment = df_result[["Barcode"] + args.bff_methods]
            classification = assignment.copy()

        valid_values = {args.singlet_str, args.negative_str, args.doublet_str}

        # Define classification function
        def classify_value(x):
            if x in valid_values:
                return x
            elif x in args.hash_list:
                return args.singlet_str
            else:
                raise ValueError(
                    f"Value '{x}' in BFF is not 'Negative', 'Doublet', or one of the hashes in the used hashes list"
                )

        classification[args.bff_methods] = classification[args.bff_methods].applymap(
            classify_value
        )

        return assignment, classification


# TODO if we keep saving AnnData/MuData in gene/hash_summary add AnnData to container for input type (https://github.com/theislab/hadge/issues/83)
# joins the assignment results with HTO, generate_anndata will return h5ad with HTO matrix
def saveAnnDataMuData(
    args: Arguments, assignment_summary: pd.DataFrame, rna_data, hto_data
):
    if args.generate_mudata or args.generate_anndata:
        assignment_summary.set_index("Barcode", inplace=True)
        hto_data.obs = hto_data.obs.join(assignment_summary, how="left").fillna(
            args.negative_str
        )

    if args.generate_anndata:
        hto_data.write(args.h5ad)

    if args.generate_mudata:
        mudata = MuData({"rna": rna_data, "hto": hto_data})
        mudata.update()
        mudata.write(args.h5mu)


def print_method_item_counts(dfs):
    """
    Takes the list of assignment/classification DataFrames (assignments/classifications) and prints a summary table:
      method name | total count | count(item1) | count(item2) | ...
    An item refers to the donor label in the assignment (HTO-1, HTO-2, ...) or the classification (singlet, doublet, negative).
    """
    rows = []
    all_items = set()

    # Extract items and their counts for every deconvolution method
    for df in dfs:
        print(df)

        method_name = df.columns[1]
        counts = df[method_name].value_counts(dropna=False)
        total = len(df)
        all_items.update(counts.index)

        row = {"method": method_name, "count_overall": total}
        row.update(counts.to_dict())
        rows.append(row)

    summary = pd.DataFrame(rows).fillna(0)

    # Convert all numeric values to int
    for col in summary.columns:
        if col != "method":
            summary[col] = summary[col].astype(int)

    # Order columns
    summary = summary[["method", "count_overall"] + sorted(list(all_items))]

    print(summary.to_string(index=False))


if __name__ == "__main__":
    # ======================== process nextflow input arguments ========================

    args = Arguments()

    # ========================= process results from modules ===========================

    assignments = []
    classifications = []

    # call all functions that process the module outputs
    functions = ProcessModuleOutput()
    function_names = list(functions.function_name_to_args_name.keys())
    for function in function_names:
        if (
            getattr(args, functions.function_name_to_args_name.get(function))
            is not None
        ):
            assignment, classification = getattr(functions, function)(args)
            assignments.append(assignment)
            classifications.append(classification)

    # ================================== save results ==================================

    # ----------------------------------- save csv's -----------------------------------

    rna_data = sc.read_10x_mtx(args.rna_matrix)
    hto_data = sc.read_10x_mtx(args.hto_matrix, gex_only=False)

    # Need to use a left join — demuxEM outputs extra barcodes not present in the input.
    # See https://github.com/lilab-bcb/demuxEM/issues/20

    # Use hto_data.obs_names() as index to perform a left join
    assignment_summary = pd.DataFrame(hto_data.obs_names, columns=["Barcode"])
    classification_summary = assignment_summary.copy()

    for assignment in assignments:
        assignment_summary = pd.merge(
            assignment_summary, assignment, on="Barcode", how="left"
        )

    for classification in classifications:
        classification_summary = pd.merge(
            classification_summary, classification, on="Barcode", how="left"
        )

    # TODO update if demuxEM works (https://github.com/theislab/hadge/issues/81)
    # .replace("", args.negative_str)
    # maybe also in demuxem()
    assignment_summary.fillna(args.negative_str).to_csv(args.assignment, index=False)
    classification_summary.fillna(args.negative_str).to_csv(
        args.classification, index=False
    )

    # -------------------------------- save mudata/anndata -----------------------------

    saveAnnDataMuData(args, assignment_summary, rna_data, hto_data)

    # -------------------------------------- versions ----------------------------------

    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "scanpy": sc.__version__,
            "numpy": np.__version__,
            "mudata": md.__version__,
            "pegasusio": io.__version__,
        }
    }

    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)
