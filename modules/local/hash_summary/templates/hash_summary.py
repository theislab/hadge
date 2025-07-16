#!/usr/bin/env python3
import pandas as pd
import scanpy as sc
# import argparse
import numpy as np
# from pathlib import Path
# from mudata import MuData
# from anndata import AnnData
# from typing import Dict
# from typing import Tuple

# parser = argparse.ArgumentParser(description="Parameters for summary process")
# parser.add_argument(
#     "--demuxem", help="Folder containing output files of demuxem", default=None
# )
# parser.add_argument(
#     "--htodemux", help="Folder containing output files of htodemux", default=None
# )
# parser.add_argument(
#     "--multiseq", help="Folder containing output files of multiseq", default=None
# )
# parser.add_argument(
#     "--hashsolo", help="Folder containing output files of hashsolo", default=None
# )
# parser.add_argument(
#     "--hashedDrops", help="Folder containing output files of hashedDrops", default=None
# )
# parser.add_argument("--bff", help="Folder containing output files of BFF", default=None)
# parser.add_argument(
#     "--gmm_demux", help="Folder containing output files of GMM-Demux", default=None
# )
# parser.add_argument("--generate_anndata", help="Generate anndata", action="store_true")
# parser.add_argument("--generate_mudata", help="Generate mudata", action="store_true")
# parser.add_argument(
#     "--read_rna_mtx",
#     help="10x-Genomics-formatted mtx directory for gene expression",
#     default=None,
# )
# parser.add_argument(
#     "--read_hto_mtx",
#     help="10x-Genomics-formatted mtx directory for HTO expression",
#     default=None,
# )
# args = parser.parse_args()


def groovy_map_str_2_dict(input_str: str) -> Dict[str, Path]:
    """
    Parses an input string in the format '[key1: value1, key2: value2, ...]' (a groovy map) into a dictionary,
    converting ALL values into Path objects.

    Args:
        input_str: A string in the specified format (e.g., '[key1: value1, key2: value2]').

    Returns:
        A dictionary where all values are Path objects.
    """
    pairs = [pair.strip() for pair in input_str.strip("[]").split(",") if pair.strip()]
    result = {}

    for pair in pairs:
        key, value = pair.split(":", 1)  # Split on first colon only
        result[key.strip()] = Path(value.strip())

    return result

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
    demuxem_res: list[str], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> None:
    assign = []
    classi = []
    params = []
    for x in demuxem_res:
        x_path = Path(x)
        obs_res_dir = find_file_with_suffix(x_path, "_obs.csv")
        obs_res = pd.read_csv(obs_res_dir)
        obs_res.rename(columns={obs_res.columns[0]: "Barcode"}, inplace=True)
        demuxem_assign = obs_res[["Barcode", "assignment"]]
        demuxem_assign.columns = ["Barcode", x_path.name]
        demuxem_assign.index = demuxem_assign.Barcode
        demuxem_assign = demuxem_assign.drop(columns=["Barcode"])
        assign.append(demuxem_assign)

        if raw_adata is not None:
            adata = raw_adata.copy()
            save_anndata(adata, demuxem_assign, x_path.name)

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            save_mudata(mudata, demuxem_assign, x_path.name)

        demuxem_classi = obs_res[["Barcode", "demux_type"]]
        demuxem_classi.columns = ["Barcode", x_path.name]
        demuxem_classi = demuxem_classi.replace("unknown", "negative")
        demuxem_classi.index = demuxem_classi.Barcode
        demuxem_classi = demuxem_classi.drop(columns=["Barcode"])
        classi.append(demuxem_classi)

        params_dir = find_file_with_suffix(x_path, "params.csv")
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [x_path.name]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("hash_summary/demuxem_assignment.csv", quoting=False)

    classi = pd.concat(classi, axis=1)
    classi.to_csv("hash_summary/demuxem_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/demuxem_params.csv")


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
    hasheddrops_res: list[str], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> None:
    assign = []
    classi = []
    params = []

    for x in hasheddrops_res:
        x_path = Path(x)
        obs_res_dir = find_file_with_suffix(x_path, "_res.csv")
        obs_res = pd.read_csv(obs_res_dir)

        obs_res["Classification"] = np.where(
            obs_res["Confident"],
            "singlet",
            np.where(obs_res["Doublet"], "doublet", "negative"),
        )
        obs_res["Best"] = np.where(
            ~obs_res["Classification"].isin(["doublet", "negative"]),
            obs_res["Best"],
            obs_res["Classification"],
        )
        obs_res.rename(columns={obs_res.columns[0]: "Barcode"}, inplace=True)

        hasheddrops_res_df = obs_res[["Barcode", "Best"]]
        hasheddrops_res_df = hasheddrops_res_df.rename(columns={"Best": x_path.name})
        assign.append(hasheddrops_res_df)

        if raw_adata is not None:
            adata = raw_adata.copy()
            save_anndata(adata, hasheddrops_res_df, x_path.name, merge_on_barcode=True)

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            save_mudata(mudata, hasheddrops_res_df, x_path.name, merge_on_barcode=True)

        hasheddrops_classi = obs_res[["Barcode", "Classification"]]
        hasheddrops_classi = hasheddrops_classi.rename(
            columns={"Classification": x_path.name}
        )
        classi.append(hasheddrops_classi)

        params_dir = find_file_with_suffix(x_path, "params.csv")
        params_res = pd.read_csv(
            params_dir, usecols=[1, 2], keep_default_na=False, index_col=0
        )
        params_res.columns = [x_path.name]
        params.append(params_res)

    assign = pd.concat(assign, axis=1).reset_index(drop=True)
    assign.to_csv("hash_summary/hasheddrops_assignment.csv", index=False, quoting=False)

    classi = pd.concat(classi, axis=1).reset_index(drop=True)
    classi.to_csv(
        "hash_summary/hasheddrops_classification.csv", index=False, quoting=False
    )

    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/hasheddrops_params.csv")


def multiseq_summary(
    results: Dict[str, Path], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> Tuple[pd.DataFrame, pd.DataFrame]:

   assignment = pd.read_csv(results['assignment'])
   assignment.columns = ["Barcode", "multiseq"]
   assignment.replace(
            {"Doublet": "doublet", "Negative": "negative"}, inplace=True
        )


        # if raw_adata is not None:
        #     print("raw_adata is not None")
        #     adata = raw_adata.copy()
        #     save_anndata(adata, multiseq_assign, x_path.name)

        # if raw_mudata is not None:
        #     mudata = raw_mudata.copy()
        #     save_mudata(mudata, multiseq_assign, x_path.name)

    classification = assignment.copy()
    classification[(classification != "doublet") & (classification != "negative")] = "singlet"

    return assignment, classification


def htodemux_summary(
    results: Dict[str, Path], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> Tuple[pd.DataFrame, pd.DataFrame]:

        assignment = pd.read_csv(results['assignment'])
        assignment.columns = ["Barcode", "htodemux"]
        assignment.replace("Doublet", "doublet", inplace=True)
        assignment.replace(
            {"Doublet": "doublet", "Negative": "negative"}, inplace=True
        )

        classification = pd.read_csv(results['classification'])
        classification.columns = ["Barcode", "htodemux"]
        classification.columns = ["Barcode", "htodemux"]
        classification.replace(
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

        return assignment, classification


def gmm_summary(
    gmmDemux_res: list[str], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> None:
    classi = []
    assign = []
    params = []
    for x in gmmDemux_res:
        x_path = Path(x)
        obs_res_dir = find_file_with_suffix(x_path, "GMM_full.csv")
        params_dir = x_path / "params.csv"
        params_res = pd.read_csv(params_dir, index_col=False)
        params_res.columns = ["Argument", x_path.name]
        params.append(params_res)

        result_row = params_res[
            params_res["Argument"].str.contains("hto_name_gmm", case=False, na=False)
        ]
        hashes_used = ""
        if not result_row.empty:
            hashes_used = result_row[x_path.name].iloc[0]
        else:
            print("No row contains the number of hashes")
        hashes = hashes_used.split(",")
        number_of_hashes = len(hashes)

        gmm_classi = pd.read_csv(obs_res_dir)
        classification_config = x_path / "GMM_full.config"
        classif_file = pd.read_csv(classification_config, header=None)

        gmm_dt = pd.DataFrame(gmm_classi)
        classification_dt = pd.DataFrame(classif_file)

        classification_dt = classification_dt.rename(
            columns={0: "Cluster_id", 1: "assignment"}
        )
        gmm_dt = gmm_dt.rename(columns={"Unnamed: 0": "Barcode"})

        def _classify_hash(row, number_hashes):
            if row == 0:
                return "negative"
            elif row > 0 and row <= number_hashes:
                print("singlet found")
                return "singlet"
            else:
                return "doublet"

        classification_dt["Classification"] = classification_dt["Cluster_id"].apply(
            lambda x: _classify_hash(x, number_of_hashes)
        )

        new_rows = []
        for _, row in gmm_dt.iterrows():
            cluster_id = row["Cluster_id"]
            matching_row_map = classification_dt[
                classification_dt["Cluster_id"] == cluster_id
            ]
            if not matching_row_map.empty:
                assignment_gmm = matching_row_map.iloc[0]["assignment"]
                classification_gmm = matching_row_map.iloc[0]["Classification"]

                new_row = {
                    "Barcode": row["Barcode"],
                    "Cluster_id": cluster_id,
                    "assignment": assignment_gmm,
                    "Classification": classification_gmm,
                }
                new_rows.append(new_row)
        merged = pd.DataFrame(new_rows)

        merged["assignment"] = merged.apply(
            lambda row: "doublet"
            if "doublet" in row["Classification"]
            else row["assignment"],
            axis=1,
        )

        gmm_dt["Classification"] = merged["Classification"]
        gmm_dt["Assignment"] = merged["assignment"]
        gmm_dt["Assignment"] = gmm_dt["Assignment"].apply(
            lambda x: "doublet" if "-" in x else x
        )
        classification_dt["Classification"] = classification_dt[
            "Classification"
        ].str.replace(" ", "")

        gmm_dt_assign = gmm_dt.drop(
            ["Cluster_id", "Confidence", "Classification"], axis=1
        )
        gmm_dt_assign["Assignment"] = gmm_dt_assign["Assignment"].str.replace(" ", "")
        gmm_dt_assign.columns = ["Barcode", x_path.name]
        assign.append(gmm_dt_assign)

        if raw_adata is not None:
            adata = raw_adata.copy()
            save_anndata(adata, gmm_dt_assign, x_path.name, merge_on_barcode=True)

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            save_mudata(mudata, gmm_dt_assign, x_path.name, merge_on_barcode=True)

        gmm_dt_classi = gmm_dt.drop(["Cluster_id", "Confidence", "Assignment"], axis=1)
        gmm_dt_classi.columns = ["Barcode", x_path.name]
        classi.append(gmm_dt_classi)

        params_dir = x_path / "params.csv"
        params_res = pd.read_csv(params_dir, index_col=False)
        params_res.columns = ["Argument", x_path.name]
        params.append(params_res)

    classi_df = pd.concat(classi, axis=1, join="outer")
    classi_df.to_csv("hash_summary/GMM_classification.csv", index=False)

    assign_df = pd.concat(assign, axis=1, join="outer")
    assign_df.to_csv("hash_summary/GMM_assignment.csv", index=False, sep=",")

    params_df = pd.concat(params, axis=1, join="outer")
    params_df.to_csv("hash_summary/GMM_params.csv", index=False)


def bff_summary(
    bff_res: list[str], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> None:
    classi = []
    assign = []
    params = []

    for x in bff_res:
        x_path = Path(x)
        obs_res_dir = find_file_with_suffix(x_path, "_bff.csv")
        bff_assign = pd.read_csv(obs_res_dir)
        data_bff = pd.DataFrame(bff_assign)
        if data_bff.empty:
            column_names = ["Barcode", x_path.name]
            df = pd.DataFrame(columns=column_names)
            classi.append(df)
            assign.append(df)
        else:
            dt_assign = data_bff.copy()
            column_names = [
                "Unnamed: 0",
                "bff_raw",
                "bff_cluster",
                "consensuscall.global",
            ]
            for column in column_names:
                if column in dt_assign.columns:
                    dt_assign = dt_assign.drop([column], axis=1)
            dt_assign.loc[dt_assign["consensuscall"] == "Doublet", "consensuscall"] = (
                "doublet"
            )
            dt_assign.loc[dt_assign["consensuscall"] == "Negative", "consensuscall"] = (
                "negative"
            )
            dt_assign["consensuscall"] = dt_assign["consensuscall"].astype("category")
            dt_assign = dt_assign.rename(
                columns={"cellbarcode": "Barcode", "consensuscall": x_path.name}
            )
            dt_assign["Barcode"] = dt_assign["Barcode"].apply(
                lambda x: x + "-1" if isinstance(x, str) else x
            )

            assign.append(dt_assign)

            if raw_adata is not None:
                adata = raw_adata.copy()
                adata.obs = adata.obs.merge(
                    dt_assign, left_index=True, right_index=True, how="left"
                )
                adata.obs.rename(columns={adata.obs.columns[0]: "donor"}, inplace=True)
                adata.obs.donor = adata.obs.donor.fillna("negative")
                adata.obs.donor = adata.obs.donor.astype(str)
                adata.write_h5ad(
                    Path("hash_summary/adata") / f"adata_with_{x_path.name}.h5ad"
                )

            if raw_mudata is not None:
                mudata = raw_mudata.copy()
                mudata["rna"].obs = mudata["rna"].obs.merge(
                    dt_assign, left_index=True, right_index=True, how="left"
                )
                mudata["rna"].obs.rename(
                    columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
                )
                mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
                mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
                mudata.update()
                mudata.write(
                    Path("hash_summary/mudata")
                    / f"mudata_with_mudata_{x_path.name}.h5mu"
                )

            dt_classi = data_bff.copy()
            column_names_class = ["bff_raw", "bff_cluster", "consensuscall"]
            for column in column_names_class:
                if column in dt_assign.columns:
                    dt_classi = dt_classi.drop([column], axis=1)
            dt_classi.loc[
                dt_classi["consensuscall.global"] == "Singlet", "consensuscall.global"
            ] = "singlet"
            dt_classi.loc[
                dt_classi["consensuscall.global"] == "Doublet", "consensuscall.global"
            ] = "doublet"
            dt_classi.loc[
                dt_classi["consensuscall.global"] == "Negative", "consensuscall.global"
            ] = "negative"
            dt_classi = dt_classi.rename(
                columns={
                    "cellbarcode": "Barcode",
                    "consensuscall.global": x_path.name,
                }
            )
            dt_classi["Barcode"] = dt_classi["Barcode"].apply(
                lambda x: x + "-1" if isinstance(x, str) else x
            )

            classi.append(dt_classi)

        params_dir = find_file_with_name(x_path, "params.csv")
        params_res = pd.read_csv(
            params_dir, usecols=[1, 2], keep_default_na=False, index_col=0
        )
        params_res.columns = [x_path.name]
        params.append(params_res)

    classi_df = pd.concat(classi, axis=1, join="outer")
    classi_df.to_csv("hash_summary/bff_classification.csv", index=False)

    assign_df = pd.concat(assign, axis=1, join="outer")
    assign_df.to_csv("hash_summary/bff_assignment.csv", index=False)

    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/bff_params.csv")


if __name__ == "__main__":
    adata = None
    mudata = None

    assignments = []
    classifications = []

    rna_data = sc.read_10x_mtx("${rna_matrix}")

    if "${generate_mudata}" == "true":
        hto_data = sc.read_10x_mtx("${hto_matrix}", gex_only=False)
        mudata = MuData({"rna": rna_data, "hto": hto_data})
        if "${generate_anndata}" == "true":
            adata = rna_data
    elif "${generate_anndata}" == "true":
        adata = rna_data

    if sum(s == "" for s in ["${htodemux_assignments}", "${htodemux_assignments}"]) == 1:
        raise ValueError("The assignment or classification file of htodemux is empty.")

    if "${htodemux_assignments}" != "":
        assignment, classification = htodemux_summary("${htodemux_assignments}", "${htodemux_classification}", adata, mudata)
        classifications.append(classification)
        assignments.append(assignment)
        #TODO use the old container again




    # if args.hashedDrops is not None:
    #     hashedDrops_res = args.hashedDrops.split(":")
    #     hasheddrops_summary(hashedDrops_res, adata, mudata)

    # if args.demuxem is not None:
    #     demuxem_res = args.demuxem.split(":")
    #     demuxem_summary(demuxem_res, adata, mudata)

    # if args.hashsolo is not None:
    #     hashsolo_res = args.hashsolo.split(":")
    #     hashsolo_summary(hashsolo_res, adata, mudata)

    # if args.multiseq is not None:
    #     multiseq_res = args.multiseq.split(":")
    #     multiseq_summary(multiseq_res, adata, mudata)



    # if args.gmm_demux is not None:
    #     gmmDemux_res = args.gmm_demux.split(":")
    #     gmm_summary(gmmDemux_res, adata, mudata)

    # if args.bff is not None:
    #     bff_res = args.bff.split(":")
    #     bff_summary(bff_res, adata, mudata)



    barcodes = rna_data.obs_names.tolist()
    assignment_summary = pd.DataFrame({'Barcodes': barcodes})
    classification_summary = pd.DataFrame({'Barcodes': barcodes})

    for assignment in assignments:
        assignment_summary = pd.merge(assignment_summary, assignment, on="Barcode", how="outer")

    assignment_summary.to_csv("${prefix}/summary_hashing_assignment.csv", index=False)

    for classification in classifications:
            classification_summary = pd.merge(classification_summary, classification, on="Barcode", how="outer")

    classification_summary.to_csv(
        "${prefix}/summary_hashing_classification.csv", index=False
    )
