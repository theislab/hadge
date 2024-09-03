#!/usr/bin/env python
import pandas as pd
import os
import scanpy as sc
import argparse
import numpy as np
from mudata import MuData

parser = argparse.ArgumentParser(description="Parameters for summary process")
parser.add_argument(
    "--demuxem", help="Folder containing output files of demuxem", default=None
)
parser.add_argument(
    "--htodemux", help="Folder containing output files of htodemux", default=None
)
parser.add_argument(
    "--multiseq", help="Folder containing output files of multiseq", default=None
)
parser.add_argument(
    "--hashsolo", help="Folder containing output files of hashsolo", default=None
)
parser.add_argument(
    "--hashedDrops", help="Folder containing output files of hashedDrops", default=None
)
parser.add_argument("--bff", help="Folder containing output files of BFF", default=None)
parser.add_argument(
    "--gmm_demux", help="Folder containing output files of GMM-Demux", default=None
)
parser.add_argument("--generate_anndata", help="Generate anndata", action="store_true")
parser.add_argument("--generate_mudata", help="Generate mudata", action="store_true")
parser.add_argument(
    "--read_rna_mtx",
    help="10x-Genomics-formatted mtx directory for gene expression",
    default=None,
)
parser.add_argument(
    "--read_hto_mtx",
    help="10x-Genomics-formatted mtx directory for HTO expression",
    default=None,
)
args = parser.parse_args()


def demuxem_summary(demuxem_res, raw_adata, raw_mudata):
    assign = []
    classi = []
    params = []
    for x in demuxem_res:
        obs_res_dir = os.path.join(
            x,
            [filename for filename in os.listdir(x) if filename.endswith("_obs.csv")][
                0
            ],
        )
        obs_res = pd.read_csv(obs_res_dir)
        obs_res.rename(columns={obs_res.columns[0]: "Barcode"}, inplace=True)
        demuxem_assign = obs_res[["Barcode", "assignment"]]
        demuxem_assign.columns = ["Barcode", os.path.basename(x)]
        # demuxem_assign.loc[:, "Barcode"] = demuxem_assign["Barcode"].apply(lambda x: x + "-1")
        # demuxem_assign["Barcode"] = demuxem_assign["Barcode"].astype(str) + "-1"
        demuxem_assign.index = demuxem_assign.Barcode
        demuxem_assign = demuxem_assign.drop(columns=["Barcode"])
        assign.append(demuxem_assign)
        if raw_adata is not None:
            print("raw_adata is not None")
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(
                demuxem_assign, left_index=True, right_index=True, how="left"
            )
            adata.obs.rename(columns={adata.obs.columns[0]: "donor"}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write(
                "hash_summary/adata/adata_with_" + os.path.basename(x) + ".h5ad"
            )

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata["rna"].obs = mudata["rna"].obs.merge(
                demuxem_assign, left_index=True, right_index=True, how="left"
            )
            mudata["rna"].obs.rename(
                columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
            )
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
            mudata.update()
            mudata.write(
                "hash_summary/mudata/mudata_with_mudata_"
                + os.path.basename(x)
                + ".h5mu"
            )

        demuxem_classi = obs_res[["Barcode", "demux_type"]]
        demuxem_classi.columns = ["Barcode", os.path.basename(x)]
        demuxem_classi = demuxem_classi.replace("unknown", "negative")
        # demuxem_classi.loc[:, "Barcode"] = demuxem_classi["Barcode"].apply(lambda x: x + '-1')
        demuxem_classi.index = demuxem_classi.Barcode
        demuxem_classi = demuxem_classi.drop(columns=["Barcode"])
        classi.append(demuxem_classi)

        params_dir = os.path.join(
            x,
            [filename for filename in os.listdir(x) if filename.endswith("params.csv")][
                0
            ],
        )
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        # params_res.rename(columns={params_res.columns[1]: os.path.basename(x)}, inplace=True)
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("hash_summary/demuxem_assignment.csv", quoting=False)

    classi = pd.concat(classi, axis=1)
    classi.to_csv("hash_summary/demuxem_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/demuxem_params.csv")


def hashsolo_summary(hashsolo_res, raw_adata, raw_mudata):
    assign = []
    classi = []
    params = []

    for x in hashsolo_res:
        obs_res_dir = os.path.join(
            x,
            [filename for filename in os.listdir(x) if filename.endswith("_res.csv")][
                0
            ],
        )
        obs_res = pd.read_csv(obs_res_dir, index_col=0)
        obs_res.index.name = "Barcode"
        hashsolo_assign = obs_res[["Classification"]]
        hashsolo_assign.columns = [os.path.basename(x)]
        hashsolo_assign = hashsolo_assign.replace(
            {"Doublet": "doublet", "Negative": "negative"}
        )
        assign.append(hashsolo_assign)

        if raw_adata is not None:
            print("raw_adata is not None")
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(
                hashsolo_assign, left_index=True, right_index=True, how="left"
            )
            adata.obs = adata.obs.rename(columns={adata.obs.columns[0]: "donor"})
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write(
                "hash_summary/adata/adata_with_" + os.path.basename(x) + ".h5ad"
            )

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata["rna"].obs = mudata["rna"].obs.merge(
                hashsolo_assign, left_index=True, right_index=True, how="left"
            )
            mudata["rna"].obs.rename(
                columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
            )
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
            mudata.update()
            mudata.write(
                "hash_summary/mudata/mudata_with_mudata_"
                + os.path.basename(x)
                + ".h5mu"
            )

        hashsolo_classi = obs_res[["most_likely_hypothesis"]]
        hashsolo_classi_copy = hashsolo_classi.copy()
        hashsolo_classi_copy["most_likely_hypothesis"] = hashsolo_classi_copy["most_likely_hypothesis"].astype(object)
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

        hashsolo_classi_copy.columns = [os.path.basename(x)]
        classi.append(hashsolo_classi_copy)

        params_dir = os.path.join(
            x,
            [filename for filename in os.listdir(x) if filename.endswith("params.csv")][
                0
            ],
        )
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("hash_summary/hashsolo_assignment.csv", quoting=False)

    classi = pd.concat(classi, axis=1)
    classi.to_csv("hash_summary/hashsolo_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/hashsolo_params.csv")


def hasheddrops_summary(hasheddrops_res, raw_adata, raw_mudata):
    assign = []
    classi = []
    params = []

    for x in hasheddrops_res:
        obs_res_dir = os.path.join(
            x,
            [filename for filename in os.listdir(x) if filename.endswith("_res.csv")][
                0
            ],
        )
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

        hasheddrops_res = obs_res[["Barcode", "Best"]]
        hasheddrops_res = hasheddrops_res.rename(columns={"Best": os.path.basename(x)})
        assign.append(hasheddrops_res)

        if raw_adata is not None:
            print("raw_adata is not None")
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(
                hasheddrops_res, left_index=True, right_on="Barcode", how="left"
            ).set_index("Barcode")
            adata.obs.rename(columns={adata.obs.columns[0]: "donor"}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write(
                "hash_summary/adata/adata_with_" + os.path.basename(x) + ".h5ad"
            )

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata["rna"].obs = (
                mudata["rna"]
                .obs.merge(
                    hasheddrops_res, left_index=True, right_on="Barcode", how="left"
                )
                .set_index("Barcode")
            )
            mudata["rna"].obs.rename(
                columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
            )
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
            mudata.update()

            mudata.write(
                "hash_summary/mudata/mudata_with_mudata_"
                + os.path.basename(x)
                + ".h5mu"
            )

        hasheddrops_classi = obs_res[["Barcode", "Classification"]]
        hasheddrops_classi = hasheddrops_classi.rename(
            columns={"Classification": os.path.basename(x)}
        )
        classi.append(hasheddrops_classi)

        params_dir = os.path.join(
            x,
            [filename for filename in os.listdir(x) if filename.endswith("params.csv")][
                0
            ],
        )
        params_res = pd.read_csv(
            params_dir, usecols=[1, 2], keep_default_na=False, index_col=0
        )
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1).reset_index(drop=True)
    assign.to_csv("hash_summary/hasheddrops_assignment.csv", index=False, quoting=False)

    classi = pd.concat(classi, axis=1).reset_index(drop=True)
    classi.to_csv(
        "hash_summary/hasheddrops_classification.csv", index=False, quoting=False
    )

    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/hasheddrops_params.csv")


def multiseq_summary(multiseq_res, raw_adata, raw_mudata):
    print("Multiseq: ",multiseq_res)
    assign = []
    params = []
    for x in multiseq_res:
        obs_res_dir = os.path.join(
            x,
            [filename for filename in os.listdir(x) if filename.endswith("_res.csv")][
                0
            ],
        )
        multiseq_assign = pd.read_csv(obs_res_dir)
        multiseq_assign.columns = ["Barcode", os.path.basename(x)]
        multiseq_assign.set_index("Barcode", inplace=True)
        multiseq_assign.replace(
            {"Doublet": "doublet", "Negative": "negative"}, inplace=True
        )

        assign.append(multiseq_assign)

        if raw_adata is not None:
            print("raw_adata is not None")
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(
                multiseq_assign, left_index=True, right_index=True, how="left"
            )
            adata.obs.rename(columns={adata.obs.columns[0]: "donor"}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write(
                "hash_summary/adata/adata_with_" + os.path.basename(x) + ".h5ad"
            )

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata["rna"].obs = mudata["rna"].obs.merge(
                multiseq_assign, left_index=True, right_index=True, how="left"
            )
            mudata["rna"].obs.rename(
                columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
            )
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
            mudata.update()
            mudata.write(
                "hash_summary/mudata/mudata_with_mudata_"
                + os.path.basename(x)
                + ".h5mu"
            )

        params_dir = os.path.join(
            x,
            [filename for filename in os.listdir(x) if filename.endswith("params.csv")][
                0
            ],
        )
        params_res = pd.read_csv(
            params_dir, usecols=[1, 2], keep_default_na=False, index_col=0
        )
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("hash_summary/multiseq_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[(classi != "doublet") & (classi != "negative")] = "singlet"
    classi.to_csv("hash_summary/multiseq_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/multiseq_params.csv")


def htodemux_summary(htodemux_res, raw_adata, raw_mudata):
    assign = []
    classi = []
    params = []
    for x in htodemux_res:
        obs_res_dir = os.path.join(
            x,
            [
                filename
                for filename in os.listdir(x)
                if filename.endswith("_assignment_htodemux.csv")
            ][0],
        )
        htodemux_assign = pd.read_csv(obs_res_dir)
        htodemux_assign.columns = ["Barcode", os.path.basename(x)]
        htodemux_assign.replace("Doublet", "doublet", inplace=True)
        htodemux_assign.replace("Negative", "negative", inplace=True)
        htodemux_assign.index = htodemux_assign.Barcode
        htodemux_assign = htodemux_assign.drop(columns=["Barcode"])
        assign.append(htodemux_assign)

        if raw_adata is not None:
            print("raw_adata is not None")
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(
                htodemux_assign, left_index=True, right_index=True, how="left"
            )
            adata.obs.rename(columns={adata.obs.columns[0]: "donor"}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write(
                "hash_summary/adata/adata_with_" + os.path.basename(x) + ".h5ad"
            )

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata["rna"].obs = (
                mudata["rna"]
                .obs.merge(
                    htodemux_assign, left_index=True, right_on="Barcode", how="left"
                )
                .set_index("Barcode")
            )
            mudata["rna"].obs.rename(
                columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
            )
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
            mudata.update()
            mudata.write(
                "hash_summary/mudata/mudata_with_mudata_"
                + os.path.basename(x)
                + ".h5mu"
            )

        obs_res_dir = os.path.join(
            x,
            [
                filename
                for filename in os.listdir(x)
                if filename.endswith("_classification_htodemux.csv")
            ][0],
        )
        htodemux_classi = pd.read_csv(obs_res_dir)
        htodemux_classi.columns = ["Barcode", os.path.basename(x)]
        htodemux_classi.replace("Singlet", "singlet", inplace=True)
        htodemux_classi.replace("Doublet", "doublet", inplace=True)
        htodemux_classi.replace("Negative", "negative", inplace=True)
        htodemux_classi.index = htodemux_classi.Barcode
        htodemux_classi = htodemux_classi.drop(columns=["Barcode"])
        classi.append(htodemux_classi)

        params_dir = os.path.join(
            x, [filename for filename in os.listdir(x) if filename == "params.csv"][0]
        )
        params_res = pd.read_csv(
            params_dir, usecols=[1, 2], keep_default_na=False, index_col=0
        )
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("hash_summary/htodemux_assignment.csv", quoting=False)

    classi = pd.concat(classi, axis=1)
    classi.to_csv("hash_summary/htodemux_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/htodemux_params.csv")


def gmm_summary(gmmDemux_res, raw_adata, raw_mudata):
    classi = []
    assign = []
    params = []
    for x in gmmDemux_res:
        obs_res_dir = os.path.join(
            x,
            [
                filename
                for filename in os.listdir(x)
                if filename.endswith("GMM_full.csv")
            ][0],
        )
        # we get the number of hashes used for the experiment from the parameters
        # so that we now how many kinds of singlets we could find in the assignment
        params_dir = os.path.join(x, "params.csv")
        params_res = pd.read_csv(params_dir, index_col=False)
        params_res.columns = ["Argument", os.path.basename(x)]
        params.append(params_res)

        ##### Get number of hashes used for the experiment
        result_row = params_res[
            params_res["Argument"].str.contains("hto_name_gmm", case=False, na=False)
        ]
        hashes_used = ""
        if not result_row.empty:
            hashes_used = result_row[os.path.basename(x)].iloc[0]
        else:
            print("No row contains the number of hashes")
        hashes = hashes_used.split(",")
        number_of_hashes = len(hashes)
        # -----------------

        # GMM assignement file - results GMM
        gmm_classi = pd.read_csv(obs_res_dir)

        # GMM full is the name given by GMM_demux per default to all results
        ##classif_file - contains the mapping of results
        classification_config = os.path.join(x, "GMM_full.config")
        classif_file = pd.read_csv(classification_config, header=None)

        # Classification and Assigment come from the same file
        gmm_dt = pd.DataFrame(gmm_classi)
        classification_dt = pd.DataFrame(classif_file)

        # change column names
        classification_dt = classification_dt.rename(
            columns={0: "Cluster_id", 1: "assignment"}
        )

        gmm_dt = gmm_dt.rename(columns={"Unnamed: 0": "Barcode"})

        # Create classification following the assignment found for the barcodes
        # we keep the original assigment and add a classification column
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

        # Compare classification guide file with classification found
        # merged = pd.merge(classification_dt, gmm_dt, on='Cluster_id', how='left')
        new_rows = []
        for index, row in gmm_dt.iterrows():
            cluster_id = row["Cluster_id"]
            matching_row_map = classification_dt[
                classification_dt["Cluster_id"] == cluster_id
            ]
            if not matching_row_map.empty:
                # Extract assignment and classification values from the matching row in df2
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
        # instead of multiple hashes, add doublet
        gmm_dt["Assignment"] = gmm_dt["Assignment"].apply(
            lambda x: "doublet" if "-" in x else x
        )
        classification_dt["Classification"] = classification_dt[
            "Classification"
        ].str.replace(" ", "")

        # Assigment for GMM-Demux
        #'Cluster_id',
        gmm_dt_assign = gmm_dt.drop(
            ["Cluster_id", "Confidence", "Classification"], axis=1
        )
        gmm_dt_assign["Assignment"] = gmm_dt_assign["Assignment"].str.replace(" ", "")
        gmm_dt_assign.columns = ["Barcode", os.path.basename(x)]
        assign.append(gmm_dt_assign)

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(
                gmm_dt_assign, left_index=True, right_on="Barcode", how="left"
            ).set_index("Barcode")
            adata.obs.rename(columns={adata.obs.columns[0]: "donor"}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write(
                "hash_summary/adata/adata_with_" + os.path.basename(x) + ".h5ad"
            )

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata["rna"].obs = (
                mudata["rna"]
                .obs.merge(
                    gmm_dt_assign, left_index=True, right_on="Barcode", how="left"
                )
                .set_index("Barcode")
            )
            mudata["rna"].obs.rename(
                columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
            )
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
            mudata.update()
            mudata.write(
                "hash_summary/mudata/mudata_with_mudata_"
                + os.path.basename(x)
                + ".h5mu"
            )

        # Classification for GMM-Demux
        gmm_dt_classi = gmm_dt.drop(["Cluster_id", "Confidence", "Assignment"], axis=1)
        gmm_dt_classi.columns = ["Barcode", os.path.basename(x)]
        classi.append(gmm_dt_classi)

        params_dir = os.path.join(x, "params.csv")
        params_res = pd.read_csv(params_dir, index_col=False)
        params_res.columns = ["Argument", os.path.basename(x)]
        params.append(params_res)

    classi_df = pd.concat(classi, axis=1, join="outer")
    classi_df.to_csv("hash_summary" + "/GMM_classification.csv", index=False)

    assign_df = pd.concat(assign, axis=1, join="outer")
    assign_df.to_csv("hash_summary" + "/GMM_assignment.csv", index=False, sep=",")

    params_df = pd.concat(params, axis=1, join="outer")
    params_df.to_csv("hash_summary" + "/GMM_params.csv", index=False)


def bff_summary(bff_res, raw_adata, raw_mudata):
    classi = []
    assign = []
    params = []

    for x in bff_res:
        obs_res_dir = os.path.join(
            x,
            [filename for filename in os.listdir(x) if filename.endswith("_bff.csv")][
                0
            ],
        )

        bff_assign = pd.read_csv(obs_res_dir)
        data_bff = pd.DataFrame(bff_assign)
        if data_bff.empty:
            # no results create empty dataframe for empty col
            column_names = ["Barcode", os.path.basename(x)]
            # Create an empty dataframe with only column names
            df = pd.DataFrame(columns=column_names)
            classi.append(df)
            assign.append(df)

        else:
            # df contain data and we save it in the same way
            dt_assign = data_bff.copy()
            # check if the columns contain results from both bff_s or only one
            column_names = [
                "Unnamed: 0",
                "bff_raw",
                "bff_cluster",
                "consensuscall.global",
            ]
            for column in column_names:
                if column in dt_assign.columns:
                    dt_assign = dt_assign.drop([column], axis=1)
            # dt_assign.loc[dt_assign["consensuscall"] == "Singlet", "consensuscall"] = "singlet"
            dt_assign.loc[
                dt_assign["consensuscall"] == "Doublet", "consensuscall"
            ] = "doublet"
            dt_assign.loc[
                dt_assign["consensuscall"] == "Negative", "consensuscall"
            ] = "negative"
            dt_assign["consensuscall"] = dt_assign["consensuscall"].astype("category")
            dt_assign = dt_assign.rename(
                columns={"cellbarcode": "Barcode", "consensuscall": os.path.basename(x)}
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
                    "hash_summary/adata/adata_with_" + os.path.basename(x) + ".h5ad"
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
                    "hash_summary/mudata/mudata_with_mudata_"
                    + os.path.basename(x)
                    + ".h5mu"
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
                    "consensuscall.global": os.path.basename(x),
                }
            )
            dt_classi["Barcode"] = dt_classi["Barcode"].apply(
                lambda x: x + "-1" if isinstance(x, str) else x
            )

            classi.append(dt_classi)

        params_dir = os.path.join(
            x, [filename for filename in os.listdir(x) if filename == "params.csv"][0]
        )
        params_res = pd.read_csv(
            params_dir, usecols=[1, 2], keep_default_na=False, index_col=0
        )
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    classi_df = pd.concat(classi, axis=1, join="outer")
    classi_df.to_csv("hash_summary" + "/bff_classification.csv", index=False)

    assign_df = pd.concat(assign, axis=1, join="outer")
    assign_df.to_csv("hash_summary" + "/bff_assignment.csv", index=False)

    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary" + "/bff_params.csv")


if __name__ == "__main__":
    adata = None
    mudata = None

    if not os.path.exists("hash_summary"):
        os.mkdir("hash_summary")

    if args.generate_anndata is True:
        os.mkdir("hash_summary/adata")
        adata = sc.read_10x_mtx(args.read_rna_mtx)

    if args.generate_mudata is True:
        os.mkdir("hash_summary/mudata")
        rna_data = sc.read_10x_mtx(args.read_rna_mtx)
        path_hto = args.read_hto_mtx
        hto_data = sc.read_10x_mtx(args.read_hto_mtx, gex_only=False)
        mudata = MuData({"rna": rna_data, "hto": hto_data})

    if args.hashedDrops is not None:
        hashedDrops_res = args.hashedDrops.split(":")
        hasheddrops_summary(hashedDrops_res, adata, mudata)

    if args.demuxem is not None:
        demuxem_res = args.demuxem.split(":")
        demuxem_summary(demuxem_res, adata, mudata)

    if args.hashsolo is not None:
        hashsolo_res = args.hashsolo.split(":")
        hashsolo_summary(hashsolo_res, adata, mudata)

    if args.multiseq is not None:
        multiseq_res = args.multiseq.split(":")
        multiseq_summary(multiseq_res, adata, mudata)

    if args.htodemux is not None:
        htodemux_res = args.htodemux.split(":")
        htodemux_summary(htodemux_res, adata, mudata)

    if args.gmm_demux is not None:
        gmmDemux_res = args.gmm_demux.split(":")
        gmm_summary(gmmDemux_res, adata, mudata)

    if args.bff is not None:
        bff_res = args.bff.split(":")
        bff_summary(bff_res, adata, mudata)

    # Read and combine assignment files
    assignment = [
        file for file in os.listdir("hash_summary") if file.endswith("_assignment.csv")
    ]
    assignment_all = pd.read_csv(os.path.join("hash_summary", assignment[0]))

    if len(assignment) > 1:
        for df in assignment[1:]:
            df = pd.read_csv(os.path.join("hash_summary", df))
            assignment_all = pd.merge(assignment_all, df, on="Barcode", how="outer")
            # print("---------assignment all--------")
            # print(assignment_all)
            # print("-----------------")

    assignment_all.to_csv("hash_summary/hashing_assignment_all.csv", index=False)

    # Read and combine classification files
    classification = [
        file
        for file in os.listdir("hash_summary")
        if file.endswith("_classification.csv")
    ]
    classification_all = pd.read_csv(os.path.join("hash_summary", classification[0]))

    if len(classification) > 1:
        for df in classification[1:]:
            df = pd.read_csv(os.path.join("hash_summary", df))
            classification_all = pd.merge(
                classification_all, df, on="Barcode", how="outer"
            )
    classification_all.to_csv(
        "hash_summary/hashing_classification_all.csv", index=False
    )
