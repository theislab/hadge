#!/usr/bin/env python
import pandas as pd
import os
import scanpy as sc
import argparse
from mudata import MuData

parser = argparse.ArgumentParser(
    description="Parameters for generating anndata and mudata"
)
parser.add_argument(
    "--assignment",
    help="Folder which contains cSV file with demultiplexing assignment",
    default=None,
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

if __name__ == "__main__":
    assignment_dir = os.path.join(
        args.assignment,
        [
            filename
            for filename in os.listdir(args.assignment)
            if filename == "all_assignment_after_match.csv"
        ][0],
    )
    if args.generate_anndata:
        if os.path.isfile(assignment_dir):
            adata = sc.read_10x_mtx(args.read_rna_mtx)
            assignment = pd.read_csv(assignment_dir, index_col=0)
            adata.obs = adata.obs.merge(
                assignment, left_index=True, right_index=True, how="left"
            )
            adata.obs = adata.obs.fillna("negative")
            adata.obs[adata.obs.columns[0]] = adata.obs[adata.obs.columns[0]].astype(
                str
            )
            adata.obs[adata.obs.columns[1]] = adata.obs[adata.obs.columns[1]].astype(
                str
            )
            adata.write("adata_with_donor_matching.h5ad")
        else:
            print(
                "Can not find assignment file, stop generating anndata. Please ensure the donor matching process has been executed to match hashing-based and genotyped-based methods."
            )

    if args.generate_mudata:
        if os.path.isfile(assignment_dir):
            rna_data = sc.read_10x_mtx(args.read_rna_mtx)
            hto_data = sc.read_10x_mtx(args.read_hto_mtx, gex_only=False)

            assignment = pd.read_csv(assignment_dir, index_col=0)
            mudata = MuData({"rna": rna_data, "hto": hto_data})

            mudata["rna"].obs = mudata["rna"].obs.merge(
                args.assignment, left_index=True, right_index=True, how="left"
            )
            mudata["rna"].obs.rename(
                columns={mudata["rna"].obs.columns[0]: "donor"}, inplace=True
            )
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.fillna("negative")
            mudata["rna"].obs.donor = mudata["rna"].obs.donor.astype(str)
            mudata.update()
            mudata.write("mudata_with_donor_matching.h5mu")
        else:
            print(
                "Can not find assignment file, stop generating anndata. Please ensure the donor matching process has been executed to match hashing-based and genotyped-based methods."
            )
