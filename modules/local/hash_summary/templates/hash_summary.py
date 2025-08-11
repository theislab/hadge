#!/usr/bin/env python3
import pandas as pd
import scanpy as sc
import numpy as np
from pathlib import Path
from mudata import MuData
from anndata import AnnData
from typing import Dict
from typing import Tuple
import pegasusio as io



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
        self.testing_inputs()

    def parse_input_args(self) -> None:


        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"

        self.rna_matrix               = "${rna_matrix}"
        self.hto_matrix               = "${hto_matrix}"
        self.htodemux_assignments     = "${htodemux_assignments}"
        self.htodemux_classification  = "${htodemux_classification}"
        self.multiseq                 = "${multiseq}"
        self.bff                      = "${bff}"
        self.demuxem                  = "${demuxem}"
        self.gmmdemux_results         = "${gmmdemux_results}"
        self.gmmdemux_config          = "${gmmdemux_config}"
        self.hasheddrops              = "${hasheddrops}"
        self.hashsolo                 = "${hashsolo}"

        self.generate_anndata         = "${generate_anndata}"
        self.generate_mudata          = "${generate_mudata}"
        self.bff_methods              = "${bff_methods}"
        self.hash_list                = "${hash_list}"

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
            "hasheddrops",
            "hashsolo",
        }

        boolean_vars = {
            "generate_anndata",
            "generate_mudata"
        }

        other_vars = {
            "bff_methods",
            "hash_list"
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
                elif input_str == "bff_methods":
                        if value_str == 'RAW':
                            return ['bff_raw']
                        elif value_str == 'CLUSTER':
                            return ['bff_cluster']
                        elif value_str == 'BOTH':
                            return ['bff_raw', 'bff_cluster','bff_consensuscall']
                        else:
                            raise ValueError(f"Methods ({value_str}) for bff not specified correctly. Choose RAW, CLUSTER or BOTH as input.")
                elif input_str == "hash_list":
                    return set(hash.strip() for hash in "${hash_list}".strip("[]").split(","))

        vars = path_vars | boolean_vars | other_vars

        for var in vars:
            raw_value = getattr(self, var)
            processed_value = _tranlate_to_python(var, raw_value)
            setattr(self, var, processed_value)

    def creat_output_dirs(self) -> None:
        directories = {
            'assignment': '_hashing_summary_assignment.csv',
            'classification': '_hashing_summary_classification.csv',
            'h5mu': '_hashing_summary.h5mu',
            'h5ad': '_hashing_summary.h5ad'
        }

        for output, directory in directories.items():
            setattr(self, output, self.prefix + directory)

    def testing_inputs(self) -> None:
        if [self.htodemux_assignments, self.htodemux_classification].count(None) == 1:
            raise ValueError("The assignment or classification file of htodemux is empty.")

        if [self.gmmdemux_results, self.gmmdemux_config].count(None) == 1:
            raise ValueError("The results or config file of gmmdemux is empty.")

    def print_args(self) -> None:
        """
        Print the arguments.
        """
        for attr in vars(self):
            print(f"{attr}: {getattr(self, attr)}")

class ProcessModuleOutput:

    def __init__(self):
        # necessary to verify which functions should to be called
        # because gmmdemux and and htodemux need two input files
        self.function_name_to_args_name = {
            'demuxem': 'demuxem',
            'hashsolo': 'hashsolo',
            'hasheddrops': 'hasheddrops',
            'multiseq': 'multiseq',
            'htodemux': 'htodemux_assignments',
            'gmmdemux': 'gmmdemux_results',
            'bff': 'bff'
        }

        self.checkHashNames = True
        self.chechEmptyInput = True
        # TODO what to do if empty assignments = []
        # T




    def demuxem(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:

        data = io.read_input(str(args.demuxem))
        classification = data.obs['demux_type'].to_frame()
        classification.reset_index(inplace=True)
        classification.columns = ["Barcode", "demuxem"]
        classification['demuxem'] = classification['demuxem'].cat.rename_categories({"unknown": args.negative_str})

        # TODO debug ob hier auch 12000
        print("debug4")
        assignment = data.obs['assignment'].to_frame()
        assignment.reset_index(inplace=True)
        assignment.columns = ["Barcode", "demuxem"]

        return assignment, classification

    def hashsolo(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:

        results = pd.read_csv(args.hashsolo, index_col=0)
        assignment = results[["Classification"]]

        assignment.columns = ["hashsolo"]
        assignment = assignment.replace(
            {"Doublet": args.doublet_str,
            "Negative": args.negative_str}
        )

        classification = results[["most_likely_hypothesis"]].copy()
        classification["most_likely_hypothesis"] = (
                classification["most_likely_hypothesis"]
                .replace({0.0: args.negative_str, 1.0: args.singlet_str, 2.0: args.doublet_str})
            )

        return assignment, classification

    def hasheddrops(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:



        # TODO remove hardcoding
        # Hardcode indexing for now
        # for later: test = pd.read_csv("hto_index_map.csv")
        idx_to_htoname_df = pd.DataFrame({
        'Index': [1, 2],
        'HTO': ['MS-11', 'MS-12']
        })

        idx_to_htoname_df.loc[len(idx_to_htoname_df)] = [np.nan, args.negative_str]
        idx_to_htoname_map = idx_to_htoname_df.set_index('Index')['HTO'].to_dict()

        obs_res = pd.read_csv(args.hasheddrops)

        obs_res["Classification"] = np.where(
            obs_res["Confident"] & obs_res["Confident"].notna(),
            args.singlet_str,
            np.where(obs_res["Doublet"] & obs_res["Doublet"].notna(), args.doublet_str, args.negative_str)
        )

        obs_res["Assignment"] = np.where(
            obs_res["Classification"].isin([args.doublet_str, args.negative_str]),
            obs_res["Classification"],
            obs_res["Best"].map(idx_to_htoname_map),
        )

        obs_res.rename(columns={obs_res.columns[0]: "Barcode"}, inplace=True)

        print(obs_res)

        classification = obs_res[["Barcode", "Classification"]].rename(columns={"Classification": "hasheddrops"})
        assignment = obs_res[["Barcode", "Assignment"]].rename(columns={"Assignment": "hasheddrops"})

        return assignment,classification

    def multiseq(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:

        assignment = pd.read_csv(args.multiseq)
        assignment.columns = ["Barcode", "multiseq"]
        assignment.replace(
                    {"Doublet": args.doublet_str, "Negative": args.negative_str}, inplace=True
                )

        classification = assignment.copy()
        classification.loc[(classification["multiseq"] != args.doublet_str) & (classification["multiseq"] != args.negative_str), "multiseq"] = args.singlet_str

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
            {"Singlet": args.singlet_str, "Doublet": args.doublet_str, "Negative": args.negative_str}, inplace=True
        )

        return assignment, classification

    def gmmdemux(self, args: Arguments) -> Tuple[pd.DataFrame, pd.DataFrame]:

        # TODO hash list is hardcoded until now
        # extend config file with Cluster_id's
        hashes = ["MS-11","MS-12"]
        number_of_hashes = len(hashes)

        df_config = pd.read_csv(args.gmmdemux_config, header=None)
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
            other=df_config["Classification"]
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

    def bff(self, args: Arguments) ->  Tuple[pd.DataFrame, pd.DataFrame]:

        # https://bimberlab.github.io/cellhashR/articles/V03-Benchmark-example.html
        #TODO used_methods =  ['RAW', 'CLUSTER', 'BOTH']
        #method = 'RAW'

        # TODO checken ob das mit den id's matched  auch für die anderen module
        # hash_list = ['ENS_ID.1','ENS_ID']

        # Load results and subset columns
        df_result = pd.read_csv(args.bff)

        df_result.rename(columns={
            'cellbarcode': 'Barcode',
            'consensuscall': 'bff_consensuscall'
        }, inplace=True)



        assignment = df_result[['Barcode'] + args.bff_methods].copy()

        # Replace 'Doublet' and 'Negative' in all used_methods columns at once
        assignment[args.bff_methods] = assignment[args.bff_methods].replace({
            'Doublet': args.doublet_str,
            'Negative': args.negative_str,
            'Discordant': 'discordant'
        })

        # Prepare sets for fast lookup
        valid_values = {args.negative_str, args.doublet_str, 'discordant'}

        # Define classification function
        def classify_value(x):
            if x in valid_values:
                return x
            elif x in args.hash_list:
                return args.singlet_str
            else:
                raise ValueError(f"Value '{x}' in BFF is not 'Negative', 'Doublet', or one of the hashes in the used hashes list")

        if len(args.bff_methods) == 3:
            # use the classification of consensuscall.global
            used_methods = args.bff_methods - ["bff_consensuscall"] + ["consensuscall.global"]
        else:
            used_methods = args.bff_methods

        # Apply classification only to used_methods columns
        classification = assignment[['Barcode'] + used_methods].copy()
        classification.rename(columns={'consensuscall.global': 'bff_consensuscall'}, inplace=True)

        classification[used_methods] = classification[used_methods].applymap(classify_value)

        return assignment, classification

def printProccedOutput() -> None:
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

    # ====================== process nextflow input arguments ======================

    args = Arguments()
    args.print_args()

    # ========================== process results from modules ==========================
    rna_data = sc.read_10x_mtx(args.rna_matrix)
    hto_data = sc.read_10x_mtx(args.hto_matrix, gex_only=False)

    # call all functions that process the module outptus and t

    assignments = []
    classifications = []

    functions = ProcessModuleOutput()
    function_names = list(functions.function_name_to_args_name.keys())
    for function in function_names:
        if getattr(args,functions.function_name_to_args_name.get(function)) is not None:
            assignment, classification = getattr(functions,function)(args)
            assignments.append(assignment)
            classifications.append(classification)





    # ========================== save results ==========================

    # ------------------------- save csv's ------------------------------

    # restructure the if statement if I keep using the hto_data
    # have to to this because demuxem has more barcodes as output that it received as input
    # https://github.com/lilab-bcb/demuxEM/issues/20

    assignment_summary = pd.DataFrame(hto_data.obs_names, columns=['Barcode'])
    classification_summary = assignment_summary.copy()

    for assignment in assignments:
        assignment_summary = pd.merge(assignment_summary, assignment, on="Barcode", how="left")

    assignment_summary.to_csv(args.assignment, index=False)

    for classification in classifications:
            classification_summary = pd.merge(classification_summary, classification, on="Barcode", how="outer")

    classification_summary.to_csv(args.classification, index=False)

    assignment_summary.set_index("Barcode", inplace=True)
    print(assignment_summary)


    # ------------------------- save mudata/anndata -----------------------
    if args.generate_mudata or args.generate_anndata:
        # join on index (Barcode)
        rna_data.obs = rna_data.obs.join(assignment_summary, how="left")
        # fill all empty of the used modules with negative values (for expression data)
        used_modules = list(assignment_summary.columns)
        for col in used_modules:
            if pd.api.types.is_categorical_dtype(rna_data.obs[col]):
                rna_data.obs[col] = rna_data.obs[col].cat.add_categories([args.negative_str])

        rna_data.obs[used_modules] = rna_data.obs[used_modules].fillna(args.negative_str)
        rna_data.obs[used_modules] = rna_data.obs[used_modules].astype(str)

        if args.generate_mudata:
            # join on index (Barcode) and create a mudata object
            hto_data.obs = hto_data.obs.join(assignment_summary, how="left")
            mudata = MuData({"rna": rna_data, "hto": hto_data})
            # mudata update?
            mudata.write(args.h5mu)

        if args.generate_anndata:
            rna_data.write(args.h5ad)
