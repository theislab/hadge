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

    def parse_input_args(self) -> None:


        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"

        self.rna_matrix               = "${rna_matrix}"
        self.hto_matrix               = "${hto_matrix}"
        self.hto_demux_assignments    = "${htodemux_assignments}"
        self.hto_demux_classification = "${htodemux_classification}"
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
            "hto_demux_assignments",
            "hto_demux_classification",
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
            # Interpret the string literal to decide if it's "[]" or something else
            if value_str.strip() == "[]":
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

    def print_args(self) -> None:
        """
        Print the arguments.
        """
        for attr in vars(self):
            print(f"{attr}: {getattr(self, attr)}")


def find_file_with_suffix(directory: Path, suffix: str) -> Path:
    return [file for file in directory.iterdir() if file.name.endswith(suffix)][0]


def find_file_with_name(directory: Path, name: str) -> Path:
    return [file for file in directory.iterdir() if file.name == name][0]


def save_anndata(
    adata: AnnData,
    assign_data: pd.DataFrame,
    basename: str,
    merge_on_barcode: bool = False,
) -> None:
    if merge_on_barcode:
        adata.obs = adata.obs.merge(
            assign_data, left_index=True, right_on="Barcode", how="left"
        ).set_index("Barcode")
    else:
        adata.obs = adata.obs.merge(
            assign_data, left_index=True, right_index=True, how="left"
        )
    adata.obs.rename(columns={adata.obs.columns[0]: "donor"}, inplace=True)
    adata.obs.donor = adata.obs.donor.fillna("negative")
    adata.obs.donor = adata.obs.donor.astype(str)
    adata.write(Path("hash_summary/adata") / f"adata_with_{basename}.h5ad")


def save_mudata(
    mudata: MuData,
    assign_data: pd.DataFrame,
    basename: str,
    merge_on_barcode: bool = False,
) -> None:
    if merge_on_barcode:
        mudata["rna"].obs = (
            mudata["rna"]
            .obs.merge(assign_data, left_index=True, right_on="Barcode", how="left")
            .set_index("Barcode")
        )
    else:
        mudata["rna"].obs = mudata["rna"].obs.merge(
            assign_data, left_index=True, right_index=True, how="left"
        )
    mudata["rna"].obs.rename(
        columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
    )
    mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
    mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
    mudata.update()
    mudata.write(Path("hash_summary/mudata") / f"mudata_with_mudata_{basename}.h5mu")


def demuxem_summary(
    results: Path, raw_adata: AnnData | None, raw_mudata: MuData | None
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    print("debug3")

    data = io.read_input(str(results))
    classi = data.obs['demux_type'].to_frame()
    classi.reset_index(inplace=True)
    classi.columns = ["Barcode", "demuxem"]
    classi['demuxem'] = classi['demuxem'].cat.rename_categories({"unknown": "negative"})

    # TODO debug ob hier auch 12000
    print("debug4")
    assign = data.obs['assignment'].to_frame()
    assign.reset_index(inplace=True)
    assign.columns = ["Barcode", "demuxem"]

    # different number of row that the other files

    # if raw_adata is not None:
    #     adata = raw_adata.copy()
    #     save_anndata(adata, demuxem_assign, x_path.name)

    # if raw_mudata is not None:
    #     mudata = raw_mudata.copy()
    #     save_mudata(mudata, demuxem_assign, x_path.name)

    return assign, classi


def hashsolo_summary(
    hashsolo_res: list[str], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> None:
    assign = []
    classi = []
    params = []

    for x in hashsolo_res:
        x_path = Path(x)
        obs_res_dir = find_file_with_suffix(x_path, "_res.csv")
        obs_res = pd.read_csv(obs_res_dir, index_col=0)
        obs_res.index.name = "Barcode"
        hashsolo_assign = obs_res[["Classification"]]
        hashsolo_assign.columns = [x_path.name]
        hashsolo_assign = hashsolo_assign.replace(
            {"Doublet": "doublet", "Negative": "negative"}
        )
        assign.append(hashsolo_assign)

        if raw_adata is not None:
            adata = raw_adata.copy()
            save_anndata(adata, hashsolo_assign, x_path.name)

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            save_mudata(mudata, hashsolo_assign, x_path.name)

        hashsolo_classi = obs_res[["most_likely_hypothesis"]]
        hashsolo_classi_copy = hashsolo_classi.copy()
        hashsolo_classi_copy["most_likely_hypothesis"] = hashsolo_classi_copy[
            "most_likely_hypothesis"
        ].astype(object)
        hashsolo_classi_copy.loc[
            hashsolo_classi_copy["most_likely_hypothesis"] == 0.0,
            "most_likely_hypothesis",
        ] = "negative"
        hashsolo_classi_copy.loc[
            hashsolo_classi_copy["most_likely_hypothesis"] == 1.0,
            "most_likely_hypothesis",
        ] = "singlet"
        hashsolo_classi_copy.loc[
            hashsolo_classi_copy["most_likely_hypothesis"] == 2.0,
            "most_likely_hypothesis",
        ] = "doublet"

        hashsolo_classi_copy.columns = [x_path.name]
        classi.append(hashsolo_classi_copy)

        params_dir = find_file_with_suffix(x_path, "params.csv")
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [x_path.name]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("hash_summary/hashsolo_assignment.csv", quoting=False)

    classi = pd.concat(classi, axis=1)
    classi.to_csv("hash_summary/hashsolo_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/hashsolo_params.csv")


def hasheddrops_summary(
    results: Path, raw_adata: AnnData | None, raw_mudata: MuData | None
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    # TODO add this to the nf-core module to have a map that indexes from the integer (in Best, etc. to the HTO name or combinations)
    # TODO if combinations is specified it the index will map to the HTO join with an +
    # otherwise it will jost take the rowname
    # # Mocking up an example dataset with 10 HTOs and 10% doublets.
    # #print(y)
    # combinations <- NULL
    # hto <- Read10X(data.dir = "/Users/luisheinzlmeier/Desktop/hto", gene.column = 2)

    # # Get the HTO names
    # hto_names <- rownames(hto)
    # if (!is.null(combinations)){
    #     hto_names <- apply(combinations, 1, function(row) paste(row, collapse = "+"))
    #     # In some applications, samples are labelled with a combination of HTOs to enable achieve greater
    #     # multiplexing throughput. This is accommodated by passing combinations to specify the valid
    #     # HTO combinations that were used for sample labelling. Each row of combinations corresponds
    #     # to a sample and should contain non-duplicated row indices of x corresponding to the HTOs used in
    #     # that sample.
    #     # Quelle: https://bioconductor.statistik.tu-dortmund.de/packages/3.18/bioc/manuals/DropletUtils/man/DropletUtils.pdf

    #     # If combinations is specified, Best instead specifies the sample (i.e., row index of combinations).
    #     # The interpretation of LogFC and LogFC2 are slightly different, and Second is not reported - see “Resolving combinatorial hashes”.
    #     # Quelle: https://rdrr.io/github/MarioniLab/DropletUtils/man/hashedDrops.html
    # }

    # # Create a data frame mapping names to indices
    # hto_map <- data.frame(
    # Index = seq_along(hto_names),
    # HTO = hto_names
    # )

    # # Write to CSV
    # write.csv(hto_map, file = "hto_index_map.csv", row.names = FALSE)

    # TODO remove hardcoding
    # Hardcode indexing for now
    # for later: test = pd.read_csv("hto_index_map.csv")
    idx_to_htoname_df = pd.DataFrame({
    'Index': [1, 2],
    'HTO': ['MS-11', 'MS-12']
    })

    idx_to_htoname_df.loc[len(idx_to_htoname_df)] = [np.nan, "negative"]
    idx_to_htoname_map = idx_to_htoname_df.set_index('Index')['HTO'].to_dict()

    obs_res = pd.read_csv(results)

    obs_res["Classification"] = np.where(
        obs_res["Confident"] & obs_res["Confident"].notna(),
        "singlet",
        np.where(obs_res["Doublet"] & obs_res["Doublet"].notna(), "doublet", "negative")
    )

    obs_res["Assignment"] = np.where(
        obs_res["Classification"].isin(["doublet", "negative"]),
        obs_res["Classification"],
        obs_res["Best"].map(idx_to_htoname_map),
    )

    obs_res.rename(columns={obs_res.columns[0]: "Barcode"}, inplace=True)

    print(obs_res)

    classi = obs_res[["Barcode", "Classification"]].rename(columns={"Classification": "hasheddrops"})
    assign = obs_res[["Barcode", "Assignment"]].rename(columns={"Assignment": "hasheddrops"})

    # if raw_adata is not None:
    #     adata = raw_adata.copy()
    #     save_anndata(adata, hasheddrops_res_df, x_path.name, merge_on_barcode=True)

    # if raw_mudata is not None:
    #     mudata = raw_mudata.copy()
    #     save_mudata(mudata, hasheddrops_res_df, x_path.name, merge_on_barcode=True)

    return assign,classi


def multiseq_summary(
    assignment: Path, raw_adata: AnnData | None, raw_mudata: MuData | None
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    assign = pd.read_csv(assignment)
    assign.columns = ["Barcode", "multiseq"]
    assign.replace(
                {"Doublet": "doublet", "Negative": "negative"}, inplace=True
            )


        # if raw_adata is not None:
        #     print("raw_adata is not None")
        #     adata = raw_adata.copy()
        #     save_anndata(adata, multiseq_assign, x_path.name)

        # if raw_mudata is not None:
        #     mudata = raw_mudata.copy()
        #     save_mudata(mudata, multiseq_assign, x_path.name)

    classi = assign.copy()
    classi.loc[(classi["multiseq"] != "doublet") & (classi["multiseq"] != "negative"), "multiseq"] = "singlet"

    return assign, classi


def htodemux_summary(
    assignment: Path, classification: Path,raw_adata: AnnData | None, raw_mudata: MuData | None
) -> Tuple[pd.DataFrame, pd.DataFrame]:

        assign = pd.read_csv(assignment)
        assign.columns = ["Barcode", "htodemux"]
        assign.replace("Doublet", "doublet", inplace=True)
        assign.replace(
            {"Doublet": "doublet", "Negative": "negative"}, inplace=True
        )

        classi = pd.read_csv(classification)
        classi.columns = ["Barcode", "htodemux"]
        classi.replace(
            {"Singlet": "singlet", "Doublet": "doublet", "Negative": "negative"}, inplace=True
        )

        # if raw_adata is not None:
        #     adata = raw_adata.copy()
        #     save_anndata(adata, assignment, "htodemux")

        # if raw_mudata is not None:
        #     mudata = raw_mudata.copy()
        #     mudata["rna"].obs = (
        #         mudata["rna"]
        #         .obs.merge(
        #             htodemux_assign, left_index=True, right_on="Barcode", how="left"
        #         )
        #         .set_index("Barcode")
        #     )
        #     mudata["rna"].obs.rename(
        #         columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
        #     )
        #     mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
        #     mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
        #     mudata.update()
        #     mudata.write(
        #         Path("hash_summary/mudata") / f"mudata_with_mudata_{x_path.name}.h5mu"
        #     )

        return assign, classi


def gmm_summary(
    results: Path, config: Path, raw_adata: AnnData | None, raw_mudata: MuData | None
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    # TODO hash list is hardcoded until now
    # extend config file with Cluster_id's
    hashes = ["MS-11","MS-12"]
    number_of_hashes = len(hashes)

    df_config = pd.read_csv(config, header=None)
    df_config.columns = ["Cluster_id", "Description"]

    def _classify_hash(cluster_id: int, number_hashes: int) -> str:
        if cluster_id == 0:
            return "negative"
        elif 1 <= cluster_id <= number_hashes:
            return "singlet"
        else:
            return "doublet"

    df_config["Classification"] = df_config["Cluster_id"].apply(
        lambda cluster_id: _classify_hash(cluster_id, number_of_hashes)
    )

    df_config["Assignment"] = df_config["Description"].where(
        df_config["Classification"] == "singlet",
        other=df_config["Classification"]
    )

    # results with Cluster_id's
    df_results = pd.read_csv(results)
    df_results.columns = ["Barcode", "Cluster_id", "Confidence"]

    df_results = df_results.merge(df_config, on="Cluster_id", how="left")

    assign = df_results[["Barcode", "Assignment"]]
    assign.columns = ["Barcode", "gmmdemux"]

    classi = df_results[["Barcode", "Classification"]]
    classi.columns = ["Barcode", "gmmdemux"]

    # if raw_adata is not None:
    #     adata = raw_adata.copy()
    #     save_anndata(adata, gmm_dt_assign, "testiii", merge_on_barcode=True)

    # if raw_mudata is not None:
    #     mudata = raw_mudata.copy()
    #     save_mudata(mudata, gmm_dt_assign, "testiii", merge_on_barcode=True)

    return assign, classi


def bff_summary(
    results: Path, method: str, hashes: set[str] ,raw_adata: AnnData | None, raw_mudata: MuData | None
) ->  Tuple[pd.DataFrame, pd.DataFrame]:

    # https://bimberlab.github.io/cellhashR/articles/V03-Benchmark-example.html
    #TODO used_methods =  ['RAW', 'CLUSTER', 'BOTH']
    #method = 'RAW'

    # TODO checken ob das mit den id's matched  auch für die anderen module
    # hash_list = ['ENS_ID.1','ENS_ID']

    if method == 'RAW':
        used_methods = ['bff_raw']
    elif method == 'CLUSTER':
        used_methods = ['bff_cluster']
    elif method == 'BOTH':
        used_methods = ['bff_raw', 'bff_cluster','bff_consensuscall']
    else:
        raise ValueError(f"Methods for bff not specified correctly. Choose RAW, CLUSTER or BOTH as input.")

    # Load results and subset columns
    df_result = pd.read_csv(results)

    df_result.rename(columns={
        'cellbarcode': 'Barcode',
        'consensuscall': 'bff_consensuscall'
    }, inplace=True)

    assign = df_result[['Barcode'] + used_methods].copy()

    # Replace 'Doublet' and 'Negative' in all used_methods columns at once
    assign[used_methods] = assign[used_methods].replace({
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
        elif x in hashes:
            return args.singlet_str
        else:
            raise ValueError(f"Value '{x}' in BFF is not 'Negative', 'Doublet', or one of the hashes in the used hashes list")

    if method == 'BOTH':
        # use the classification of consensuscall.global
        used_methods = used_methods - ["bff_consensuscall"] + ["consensuscall.global"]

    # Apply classification only to used_methods columns
    classi = assign[['Barcode'] + used_methods].copy()
    classi.rename(columns={'consensuscall.global': 'bff_consensuscall'}, inplace=True)

    classi[used_methods] = classi[used_methods].applymap(classify_value)

    return assign, classi

            # if raw_adata is not None:
        #     adata = raw_adata.copy()
        #     adata.obs = adata.obs.merge(
        #         dt_assign, left_index=True, right_index=True, how="left"
        #     )
        #     adata.obs.rename(columns={adata.obs.columns[0]: "donor"}, inplace=True)
        #     adata.obs.donor = adata.obs.donor.fillna("negative")
        #     adata.obs.donor = adata.obs.donor.astype(str)
        #     adata.write_h5ad(
        #         Path("hash_summary/adata") / f"adata_with_{x_path.name}.h5ad"
        #     )

        # if raw_mudata is not None:
        #     mudata = raw_mudata.copy()
        #     mudata["rna"].obs = mudata["rna"].obs.merge(
        #         dt_assign, left_index=True, right_index=True, how="left"
        #     )
        #     mudata["rna"].obs.rename(
        #         columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
        #     )
        #     mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
        #     mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
        #     mudata.update()
        #     mudata.write(
        #         Path("hash_summary/mudata")
        #         / f"mudata_with_mudata_{x_path.name}.h5mu"
        #     )


if __name__ == "__main__":

    # parse and print arguments
    args = Arguments()
    args.print_args()

    print("debug1")

    adata = None
    mudata = None

    assignments = []
    classifications = []

    #hashes = set(ast.literal_eval("${hash_list}"))

    hashes = set(hash.strip() for hash in "${hash_list}".strip("[]").split(","))

    rna_data = sc.read_10x_mtx("${rna_matrix}")
    hto_data = sc.read_10x_mtx("${hto_matrix}", gex_only=False)

    if sum(s == "" for s in ["${htodemux_assignments}", "${htodemux_assignments}"]) == 1:
        raise ValueError("The assignment or classification file of htodemux is empty.")

    if "${htodemux_assignments}" != "":
        assignment, classification = htodemux_summary(Path("${htodemux_assignments}"), Path("${htodemux_classification}"), adata, mudata)
        assignments.append(assignment)
        classifications.append(classification)
        #TODO use the old container again

    if "${multiseq}" != "":
        assignment, classification = multiseq_summary(Path("${multiseq}"), adata, mudata)
        assignments.append(assignment)
        classifications.append(classification)

    if "${demuxem}" != "":
        assignment, classification = demuxem_summary(Path("${demuxem}"), adata, mudata)
        assignments.append(assignment)
        classifications.append(classification)

    if "${hasheddrops}" != "":
        assignment, classification = hasheddrops_summary(Path("${hasheddrops}"), adata, mudata)
        assignments.append(assignment)
        classifications.append(classification)

    if sum(s == "" for s in ["${gmmdemux_results}", "${gmmdemux_config}"]) == 1:
        raise ValueError("The assignment or classification file of htodemux is empty.")

    if "${gmmdemux_results}" != "":
        assignment, classification = gmm_summary(Path("${gmmdemux_results}"), Path("${gmmdemux_config}"), adata, mudata)
        assignments.append(assignment)
        classifications.append(classification)

    if "${bff}" != "":
        assignment, classification = bff_summary(Path("${bff}"), "${bff_methods}", hashes, adata, mudata)
        assignments.append(assignment)
        classifications.append(classification)

    if "${hashsolo}" != "":
        assignment, classification = hashsolo_summary(Path("${hashsolo}"), adata, mudata)
        assignments.append(assignment)
        classifications.append(classification)

    # TODO what to do if empty assignments = []

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

    # assignment_summary = assignments.pop(0)
    # classification_summary = classifications.pop(0)

    # restructure the if statement if I keep using the hto_data
    # have to to this because demuxem has more barcodes as output that it received as input
    # https://github.com/lilab-bcb/demuxEM/issues/20
    hto_data = sc.read_10x_mtx("${hto_matrix}", gex_only=False)
    assignment_summary = pd.DataFrame(hto_data.obs_names, columns=['Barcode'])
    classification_summary = assignment_summary.copy()

    for assignment in assignments:
        assignment_summary = pd.merge(assignment_summary, assignment, on="Barcode", how="left")

    assignment_summary.to_csv("${prefix}_hashing_summary_assignment.csv", index=False)

    for classification in classifications:
            classification_summary = pd.merge(classification_summary, classification, on="Barcode", how="outer")

    classification_summary.to_csv(
        "${prefix}_hashing_summary_classification.csv", index=False
    )


    assignment_summary = assignment_summary.set_index("Barcode", inplace=True)
    print(assignment_summary)

    if "${generate_mudata}" == "true" or "${generate_anndata}" == "true":
        # join on index (Barcode)
        rna_data.obs = rna_data.obs.join(assignment_summary, how="left")
        # fill all empty of the used modules with negative values (for expression data)
        used_modules = assignment_summary.column_names
        rna_data.obs[used_modules] = rna_data.obs[used_modules].fillna("negative")
        rna_data.obs[used_modules] = rna_data.obs[used_modules].astype(str)

        if "${generate_mudata}" == "true":
            # join on index (Barcode) and create a mudata object
            hto_data.obs = hto_data.obs.join(assignment_summary, how="left")
            mudata = MuData({"rna": rna_data, "hto": hto_data})
            mudata.write("${prefix}_hashing_summary.h5mu")

        if "${generate_anndata}" == "true":
            rna_data.write("${prefix}_hashing_summary.h5ad")
