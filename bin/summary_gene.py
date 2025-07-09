#!/usr/bin/env python
import argparse
import numpy as np
import scanpy as sc
import pandas as pd
from pathlib import Path
from mudata import MuData
from anndata import AnnData


parser = argparse.ArgumentParser(description="Parameters for summary process")
parser.add_argument(
    "--demuxlet", help="Folder containing output files of Demuxlet", default=None
)
parser.add_argument(
    "--freemuxlet", help="Folder containing output files of Freemuxlet", default=None
)
parser.add_argument(
    "--vireo", help="Folder containing output files of Vireo", default=None
)
parser.add_argument(
    "--souporcell", help="Folder containing output files of Souporcell", default=None
)
parser.add_argument(
    "--scsplit", help="Folder containing output files of scSplit", default=None
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


def find_file_with_suffix(directory: Path, suffix: str) -> Path:
    return [file for file in directory.iterdir() if file.name.endswith(suffix)][0]


def find_file_with_name(directory: Path, name: str) -> Path:
    for root, dirs, files in directory.rglob("*"):
        if name in files:
            return root / name
    return Path("")


def process_assignment(obs_res: pd.DataFrame, basename: str) -> pd.DataFrame:
    obs_res.rename(columns={"BARCODE": "Barcode", "Assignment": basename}, inplace=True)
    obs_res.set_index("Barcode", inplace=True)
    return obs_res[[basename]]


def save_anndata(adata: AnnData, assign_data: pd.DataFrame, basename: str) -> None:
    adata.obs = adata.obs.merge(
        assign_data, left_index=True, right_index=True, how="left"
    )
    adata.obs.rename(columns={adata.obs.columns[0]: "donor"}, inplace=True)
    adata.obs.donor = adata.obs.donor.fillna("negative")
    adata.obs.donor = adata.obs.donor.astype(str)
    adata.write(Path("genetic_summary/adata") / f"adata_with_{basename}.h5ad")


def save_mudata(mudata: MuData, assign_data: pd.DataFrame, basename: str) -> None:
    mudata["rna"].obs = (
        mudata["rna"]
        .obs.merge(assign_data, left_index=True, right_on="Barcode", how="left")
        .set_index("Barcode")
    )
    mudata["rna"].obs.rename(
        columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
    )
    mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
    mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
    mudata.update()
    mudata.write(Path("genetic_summary/mudata") / f"mudata_with_{basename}.h5mu")


def demuxlet_summary(
    demuxlet_res: list[str], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> None:
    assign = []
    params = []

    for res in demuxlet_res:
        x_path = Path(res)
        obs_res_dir = find_file_with_suffix(x_path, ".best")
        obs_res = pd.read_csv(obs_res_dir, sep="\t")
        obs_res = obs_res.iloc[:, [1, 4, 5]]
        obs_res["Assignment"] = np.where(
            obs_res["BEST.GUESS"].str.split(",").str[0]
            == obs_res["BEST.GUESS"].str.split(",").str[1],
            obs_res["BEST.GUESS"].str.split(",").str[0],
            "doublet",
        )
        obs_res["Assignment"] = np.where(
            obs_res["DROPLET.TYPE"] == "AMB", "negative", obs_res["Assignment"]
        )
        demuxlet_assign = process_assignment(obs_res, x_path.name)

        if raw_adata is not None:
            adata = raw_adata.copy()
            save_anndata(adata, demuxlet_assign, x_path.name)
        assign.append(demuxlet_assign)

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            save_mudata(mudata, demuxlet_assign, x_path.name)

        params_dir = find_file_with_suffix(x_path, "params.csv")
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [x_path.name]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("genetic_summary/demuxlet_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[~classi.isin(["doublet", "negative"])] = "singlet"
    classi.to_csv("genetic_summary/demuxlet_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("genetic_summary/demuxlet_params.csv")


def freemuxlet_summary(
    freemuxlet_res: list[str], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> None:
    assign = []
    params = []

    for res in freemuxlet_res:
        x_path = Path(res)
        obs_res_dir = find_file_with_suffix(x_path, ".clust1.samples.gz")
        obs_res = pd.read_csv(obs_res_dir, sep="\t")
        obs_res = obs_res.iloc[:, [1, 4, 5]]
        obs_res["Assignment"] = np.where(
            obs_res["BEST.GUESS"].str.split(",").str[0]
            == obs_res["BEST.GUESS"].str.split(",").str[1],
            obs_res["BEST.GUESS"].str.split(",").str[0],
            "doublet",
        )
        obs_res["Assignment"] = np.where(
            obs_res["DROPLET.TYPE"] == "AMB", "negative", obs_res["Assignment"]
        )
        freemuxlet_assign = process_assignment(obs_res, x_path.name)

        if raw_adata is not None:
            adata = raw_adata.copy()
            save_anndata(adata, freemuxlet_assign, x_path.name)

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            save_mudata(mudata, freemuxlet_assign, x_path.name)

        assign.append(freemuxlet_assign)

        params_dir = find_file_with_suffix(x_path, "params.csv")
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [x_path.name]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("genetic_summary/freemuxlet_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[~classi.isin(["doublet", "negative"])] = "singlet"
    classi.to_csv("genetic_summary/freemuxlet_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("genetic_summary/freemuxlet_params.csv")


def souporcell_summary(
    souporcell_res: list[str], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> None:
    assign = []
    params = []

    for res in souporcell_res:
        x_path = Path(res)
        obs_res_dir = find_file_with_name(x_path, "clusters.tsv")
        obs_res = pd.read_csv(obs_res_dir, sep="\t")
        obs_res = obs_res.iloc[:, 0:3]
        obs_res.loc[obs_res["status"] == "doublet", "assignment"] = "doublet"
        obs_res.loc[obs_res["status"] == "unassigned", "assignment"] = "negative"
        obs_res.rename(
            columns={"barcode": "Barcode", "assignment": x_path.name},
            inplace=True,
        )
        obs_res.set_index("Barcode", inplace=True)
        obs_res = obs_res[[x_path.name]]

        if raw_adata is not None:
            adata = raw_adata.copy()
            save_anndata(adata, obs_res, x_path.name)

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            save_mudata(mudata, obs_res, x_path.name)

        assign.append(obs_res)

        params_dir = find_file_with_suffix(x_path, "params.csv")
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [x_path.name]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("genetic_summary/souporcell_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[~classi.isin(["doublet", "negative"])] = "singlet"
    classi.to_csv("genetic_summary/souporcell_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("genetic_summary/souporcell_params.csv")


def vireo_summary(
    vireo_res: list[str], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> None:
    assign = []
    params = []

    for res in vireo_res:
        x_path = Path(res)
        obs_res_dir = find_file_with_name(x_path, "donor_ids.tsv")
        obs_res = pd.read_csv(obs_res_dir, sep="\t")
        obs_res.iloc[:, [0, 1]]
        obs_res[obs_res == "unassigned"] = "negative"
        obs_res.rename(
            columns={"cell": "Barcode", "donor_id": x_path.name}, inplace=True
        )
        obs_res.set_index("Barcode", inplace=True)
        obs_res = obs_res[[x_path.name]]

        if raw_adata is not None:
            adata = raw_adata.copy()
            save_anndata(adata, obs_res, x_path.name)

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            save_mudata(mudata, obs_res, x_path.name)

        assign.append(obs_res)

        params_dir = find_file_with_suffix(x_path, "params.csv")
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [x_path.name]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("genetic_summary/vireo_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[~classi.isin(["doublet", "negative"])] = "singlet"
    classi.to_csv("genetic_summary/vireo_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("genetic_summary/vireo_params.csv")


def scsplit_summary(
    scsplit_res: list[str], raw_adata: AnnData | None, raw_mudata: MuData | None
) -> None:
    assign = []
    params = []

    for x in scsplit_res:
        x_path = Path(x)
        obs_res_dir = find_file_with_name(x_path, "scSplit_result.csv")
        obs_res = pd.read_table(obs_res_dir)
        obs_res["Assignment"] = obs_res["Cluster"].str.split("-").str[1]
        obs_res["Classification"] = obs_res["Cluster"].str.split("-").str[0]
        obs_res.loc[obs_res["Classification"] == "DBL", "Assignment"] = "doublet"
        obs_res = obs_res.drop(columns=["Cluster", "Classification"])
        obs_res.set_index("Barcode", inplace=True)
        obs_res.columns = [x_path.name]

        if raw_adata is not None:
            adata = raw_adata.copy()
            save_anndata(adata, obs_res, x_path.name)

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            save_mudata(mudata, obs_res, x_path.name)

        assign.append(obs_res)

        params_dir = find_file_with_suffix(x_path, "params.csv")
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [x_path.name]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("genetic_summary/scsplit_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[(classi != "negative") & (classi != "doublet")] = "singlet"
    classi.to_csv("genetic_summary/scsplit_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("genetic_summary/scsplit_params.csv")


if __name__ == "__main__":
    adata = None
    mudata = None
    genetic_summary = Path("genetic_summary")
    genetic_summary.mkdir(exist_ok=True)

    if args.generate_anndata is True:
        (genetic_summary / "adata").mkdir(exist_ok=True)
        adata = sc.read_10x_mtx(args.read_rna_mtx)

    if args.generate_mudata is True:
        (genetic_summary / "mudata").mkdir(exist_ok=True)
        rna_data = sc.read_10x_mtx(args.read_rna_mtx)
        hto_data = sc.read_10x_mtx(args.read_hto_mtx, gex_only=False)
        mudata = MuData({"rna": rna_data, "hto": hto_data})

    if args.demuxlet is not None:
        demuxlet_res = args.demuxlet.split(":")
        demuxlet_summary(demuxlet_res, adata, mudata)
        print("Demuxlet result found")

    if args.freemuxlet is not None:
        freemuxlet_res = args.freemuxlet.split(":")
        freemuxlet_summary(freemuxlet_res, adata, mudata)
        print("Freemuxlet result found")

    if args.vireo is not None:
        vireo_res = args.vireo.split(":")
        vireo_summary(vireo_res, adata, mudata)
        print("Vireo result found")

    if args.scsplit is not None:
        scsplit_res = args.scsplit.split(":")
        scsplit_summary(scsplit_res, adata, mudata)
        print("scSplit result found")

    if args.souporcell is not None:
        souporcell_res = args.souporcell.split(":")
        souporcell_summary(souporcell_res, adata, mudata)
        print("Souporcell result found")

    assignment = [
        file
        for file in genetic_summary.iterdir()
        if file.name.endswith("_assignment.csv")
    ]
    assignment_all = pd.read_csv(genetic_summary / assignment[0])

    if len(assignment) > 1:
        for df in assignment[1:]:
            df = pd.read_csv(genetic_summary / df)
            assignment_all = pd.merge(assignment_all, df, on="Barcode", how="outer")
    assignment_all.to_csv("genetic_summary/genetic_assignment_all.csv", index=False)

    classification = [
        file
        for file in genetic_summary.iterdir()
        if file.name.endswith("_classification.csv")
    ]
    classification_all = pd.read_csv(genetic_summary / classification[0])

    if len(classification) > 1:
        for df in classification[1:]:
            df = pd.read_csv(genetic_summary / df)
            classification_all = pd.merge(
                classification_all, df, on="Barcode", how="outer"
            )
    classification_all.to_csv(
        "genetic_summary/genetic_classification_all.csv", index=False
    )
