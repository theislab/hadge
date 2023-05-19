#!/usr/bin/env python
import pandas as pd
import os
import scanpy as sc
import argparse
import muon as mu

parser = argparse.ArgumentParser(description="Parameters for generating anndata and mudata")
parser.add_argument("--assignment", help="Folder which contains cSV file with demultiplexing assignment", default=None)
parser.add_argument("--generate_anndata", help="Generate anndata", action='store_true')
parser.add_argument("--generate_mudata", help="Generate mudata", action='store_true')
parser.add_argument("--read_rna_mtx", help="10x-Genomics-formatted mtx directory for gene expression", default=None)
parser.add_argument("--read_hto_mtx", help="10x-Genomics-formatted mtx directory for HTO expression", default=None)

args = parser.parse_args()

if __name__ == '__main__':
    if args.generate_anndata:
        adata = sc.read_10x_mtx(args.read_rna_mtx)
        assignment_dir = os.path.join(args.assignment, 
                                    [filename for filename in os.listdir(args.assignment) if filename == "all_assignment_after_match.csv"][0])

        assignment = pd.read_csv(assignment_dir, index_col = 0)
        adata.obs = adata.obs.merge(assignment, left_index=True, right_index=True, how='left')
        adata.obs= adata.obs.fillna("negative")
        adata.obs[adata.obs.columns[0]]= adata.obs[adata.obs.columns[0]].astype(str)
        adata.obs[adata.obs.columns[1]]= adata.obs[adata.obs.columns[1]].astype(str)
        adata.write("adata_with_donor_matching.h5ad")
    
    if args.generate_mudata:
        # write mudata_with_donor_matching.h5mu data
        pass
