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
import pegasusio as io

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
        self.testing_inputs()

    def parse_input_args(self) -> None:
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"

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

        self.bff_methods = "${bff_methods}"
        self.hash_list = "${hash_list}"

        path_vars = {
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

        other_vars = {"bff_methods", "hash_list"}

        def _tranlate_to_python(input_str, value_str):
            if value_str.strip() == "":
                return None
            else:
                if input_str in path_vars:
                    path = Path(value_str)
                    if not path.exists():
                        raise FileNotFoundError(f"Path does not exist: {path}")
                    return path
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

        vars = path_vars | other_vars

        for var in vars:
            raw_value = getattr(self, var)
            processed_value = _tranlate_to_python(var, raw_value)
            setattr(self, var, processed_value)

    def creat_output_dirs(self) -> None:
        directories = {
            "assignment": "_hashing_summary_assignment.csv",
            "classification": "_hashing_summary_classification.csv",
            "overview_assignment": "_hashing_overview_assignment.csv",
            "overview_classification": "_hashing_overview_classification.csv",
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
        df = data.obs[["assignment","demux_type"]].copy()
        df.index.name = "Barcode"
        df.reset_index(inplace=True)
        df["demux_type"] = df["demux_type"].cat.rename_categories(
            lambda x: args.negative_str if x == "unknown" else x
        )
        df["assignment"] = df["assignment"].cat.add_categories([args.negative_str, args.doublet_str])
        df.loc[df["demux_type"] == args.negative_str, "assignment"] = args.negative_str
        df.loc[df["demux_type"] == args.doublet_str, "assignment"] = args.doublet_str
        df["assignment"] = df["assignment"].cat.remove_unused_categories()
        assignment = df[["Barcode", "assignment"]].rename(columns={"assignment": "demuxem"})
        classification = df[["Barcode", "demux_type"]].rename(columns={"demux_type": "demuxem"})
        # TODO demuxem: demuxem has more output barcodes than input barcodes metioned here: https://github.com/lilab-bcb/demuxEM/issues/20
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

    # save overview tables
    overview_assignment = create_overview_table(assignments)
    overview_assignment.to_csv(args.overview_assignment, index=False)
    overview_classifications = create_overview_table(classifications)
    overview_classifications.to_csv(args.overview_classification, index=False)

    # save summary of all deconvolution methods
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

    for df in [assignment_summary, classification_summary]:
        for col in df.select_dtypes(["category"]):
            if args.negative_str not in df[col].cat.categories:
                df[col] = df[col].cat.add_categories(args.negative_str)

    # TODO demuxem: update if demuxEM works (https://github.com/theislab/hadge/issues/81)
    # .replace("", args.negative_str)
    # maybe also in demuxem()
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
            "pegasusio": io.__version__,
        }
    }

    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)
