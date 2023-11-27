#!/usr/bin/env python
import os
import argparse
import numpy as np
import scanpy as sc
import pandas as pd
from mudata import MuData


parser = argparse.ArgumentParser(description="Parameters for summary process")
parser.add_argument("--demuxlet", help="Folder containing output files of Demuxlet", default=None)
parser.add_argument("--freemuxlet", help="Folder containing output files of Freemuxlet", default=None)
parser.add_argument("--vireo", help="Folder containing output files of Vireo", default=None)
parser.add_argument("--souporcell", help="Folder containing output files of Souporcell", default=None)
parser.add_argument("--scsplit", help="Folder containing output files of scSplit", default=None)
parser.add_argument("--generate_anndata", help="Generate anndata", action='store_true')
parser.add_argument("--generate_mudata", help="Generate mudata", action='store_true')
parser.add_argument("--read_rna_mtx", help="10x-Genomics-formatted mtx directory for gene expression", default=None)
parser.add_argument("--read_hto_mtx", help="10x-Genomics-formatted mtx directory for HTO expression", default=None)

args = parser.parse_args()

def demuxlet_summary(demuxlet_res, raw_adata, raw_mudata):
    assign = []
    params = []
    for x in demuxlet_res:
        obs_res_dir = [file for file in os.listdir(x) if file.endswith('.best')][0]
        obs_res = pd.read_csv(os.path.join(x, obs_res_dir), sep='\t')
        obs_res = obs_res.iloc[:, [1, 4, 5]]
        obs_res['Assignment'] = np.where(obs_res['BEST.GUESS'].str.split(',').str[0] == obs_res['BEST.GUESS'].str.split(',').str[1], 
                                         obs_res['BEST.GUESS'].str.split(',').str[0], "doublet")
        obs_res['Assignment'] = np.where(obs_res['DROPLET.TYPE'] == 'AMB', 'negative', obs_res['Assignment'])
        obs_res.rename(columns={"BARCODE": "Barcode", "Assignment": os.path.basename(x)}, inplace=True)
        obs_res.set_index('Barcode', inplace=True)
        demuxlet_assign = obs_res[[os.path.basename(x)]]
        
        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(demuxlet_assign, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("genetic_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")
        assign.append(demuxlet_assign)

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata['rna'].obs = mudata['rna'].obs.merge(demuxlet_assign, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
            mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
            mudata.update()
            mudata.write("genetic_summary/mudata/mudata_with_"+ os.path.basename(x)+".h5mu") 

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("genetic_summary/demuxlet_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[~classi.isin(["doublet", "negative"])] = "singlet"
    classi.to_csv("genetic_summary/demuxlet_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("genetic_summary/demuxlet_params.csv")

def freemuxlet_summary(freemuxlet_res, raw_adata, raw_mudata):
    assign = []
    params = []

    for x in freemuxlet_res:
        obs_res_dir = [file for file in os.listdir(x) if file.endswith('.clust1.samples.gz')][0]
        obs_res = pd.read_csv(os.path.join(x, obs_res_dir), sep='\t')
        obs_res = obs_res.iloc[:, [1, 4, 5]]
        obs_res['Assignment'] = np.where(obs_res['BEST.GUESS'].str.split(',').str[0] == obs_res['BEST.GUESS'].str.split(',').str[1], 
                                                obs_res['BEST.GUESS'].str.split(',').str[0], "doublet")
        obs_res['Assignment'] = np.where(obs_res['DROPLET.TYPE'] == 'AMB', 'negative', obs_res['Assignment'])
        obs_res.rename(columns={"BARCODE": "Barcode", "Assignment": os.path.basename(x)}, inplace=True)
        obs_res.set_index('Barcode', inplace=True)
        freemuxlet_assign = obs_res[[os.path.basename(x)]]
        
        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(freemuxlet_assign, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("genetic_summary/adata/adata_with_"+ os.path.basename(x)+".h5ad")

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata['rna'].obs = mudata['rna'].obs.merge(freemuxlet_assign, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
            mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
            mudata.update()
            mudata.write("genetic_summary/mudata/mudata_with_"+ os.path.basename(x)+".h5mu")

        assign.append(freemuxlet_assign)

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("genetic_summary/freemuxlet_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[~classi.isin(["doublet", "negative"])] = "singlet"
    classi.to_csv("genetic_summary/freemuxlet_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("genetic_summary/freemuxlet_params.csv")

def souporcell_summary(souporcell_res, raw_adata, raw_mudata):
    assign = []
    params = []
    for x in souporcell_res:
        obs_res_dir = ""
        for root, dirs, files in os.walk(x):
            if "clusters.tsv" in files:
                obs_res_dir = os.path.join(root, "clusters.tsv")
        obs_res = pd.read_csv(os.path.join(x, obs_res_dir), sep='\t')
        obs_res = obs_res.iloc[:, 0:3]
        obs_res.loc[obs_res['status'] == 'doublet', 'assignment'] = 'doublet'
        obs_res.loc[obs_res['status'] == 'unassigned', 'assignment'] = 'negative'
        obs_res.rename(columns={'barcode': 'Barcode', 'assignment': os.path.basename(x)}, inplace=True)
        obs_res.set_index('Barcode', inplace=True)
        obs_res = obs_res[[os.path.basename(x)]]

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(obs_res, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("genetic_summary/adata/adata_with_"+ os.path.basename(x)+".h5ad")

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata['rna'].obs = mudata['rna'].obs.merge(obs_res, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
            mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
            mudata.update()
            mudata.write("genetic_summary/mudata/mudata_with_"+ os.path.basename(x)+".h5mu")

        assign.append(obs_res)

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("genetic_summary/souporcell_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[~classi.isin(["doublet", "negative"])] = "singlet"
    classi.to_csv("genetic_summary/souporcell_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("genetic_summary/souporcell_params.csv")

def vireo_summary(vireo_res, raw_adata, raw_mudata):
    assign = []
    params = []

    for x in vireo_res:
        obs_res_dir = ""
        for root, dirs, files in os.walk(x):
            if "donor_ids.tsv" in files:
                obs_res_dir = os.path.join(root, "donor_ids.tsv")
        obs_res = pd.read_csv(os.path.join(x, obs_res_dir), sep='\t')
        bs_res = obs_res.iloc[:, [0, 1]]
        obs_res[obs_res == "unassigned"] = "negative"
        obs_res.rename(columns={'cell': 'Barcode', 'donor_id': os.path.basename(x)}, inplace=True)
        obs_res.set_index('Barcode', inplace=True)
        obs_res = obs_res[[os.path.basename(x)]]

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(obs_res, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("genetic_summary/adata/adata_with_"+ os.path.basename(x)+".h5ad")

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata['rna'].obs = mudata['rna'].obs.merge(obs_res, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
            mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
            mudata.update()
            mudata.write("genetic_summary/mudata/mudata_with_"+ os.path.basename(x)+".h5mu")

        assign.append(obs_res)

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("genetic_summary/vireo_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[~classi.isin(["doublet", "negative"])] = "singlet"
    classi.to_csv("genetic_summary/vireo_classification.csv", quoting=False)

    params = pd.concat(params, axis=1)
    params.to_csv("genetic_summary/vireo_params.csv")

def scsplit_summary(scsplit_res, raw_adata, raw_mudata):
    assign = []
    params = []

    for x in scsplit_res:
        obs_res_dir = next((os.path.join(root, "scSplit_result.csv") for root, dirs, files in os.walk(x) if "scSplit_result.csv" in files),"")
        obs_res = pd.read_table(obs_res_dir)
        obs_res['Assignment'] = obs_res['Cluster'].str.split('-').str[1]
        obs_res['Classification'] = obs_res['Cluster'].str.split('-').str[0]
        obs_res.loc[obs_res['Classification'] == 'DBL', 'Assignment'] = 'doublet'        
        obs_res = obs_res.drop(columns=['Cluster', 'Classification'])
        obs_res.set_index('Barcode', inplace=True)
        obs_res.columns = [os.path.basename(x)]
        
        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(obs_res, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("genetic_summary/adata/adata_with_"+ os.path.basename(x)+".h5ad")
        
        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata['rna'].obs = mudata['rna'].obs.merge(obs_res, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
            mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
            mudata.update()
            mudata.write("genetic_summary/mudata/mudata_with_"+ os.path.basename(x)+".h5mu")

        assign.append(obs_res)

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("genetic_summary/scsplit_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[(classi != 'negative') & (classi != 'doublet')] = 'singlet'
    classi.to_csv("genetic_summary/scsplit_classification.csv", quoting=False)
   
    params = pd.concat(params, axis=1)
    params.to_csv("genetic_summary/scsplit_params.csv")

if __name__ == '__main__':
    adata = None
    mudata = None
    if not os.path.exists("genetic_summary"):
        os.mkdir("genetic_summary")
    
    if args.generate_anndata is True:
        os.mkdir("genetic_summary/adata")
        adata = sc.read_10x_mtx(args.read_rna_mtx)

    if args.generate_mudata is True:
        if not os.path.exists("genetic_summary/mudata"):
            os.mkdir("genetic_summary/mudata")
        rna_data = sc.read_10x_mtx(args.read_rna_mtx)
        hto_data = sc.read_10x_mtx(args.read_hto_mtx, gex_only=False)
        mudata = MuData({"rna": rna_data, "hto": hto_data })

    if args.demuxlet is not None:
        demuxlet_res = args.demuxlet.split(':')
        demuxlet_summary(demuxlet_res, adata, mudata)
        print("Demuxlet result found")

    if args.freemuxlet is not None:
        freemuxlet_res = args.freemuxlet.split(':')
        freemuxlet_summary(freemuxlet_res, adata, mudata)
        print("Freemuxlet result found")

    if args.vireo is not None:
        vireo_res = args.vireo.split(':')
        vireo_summary(vireo_res, adata, mudata)
        print("Vireo result found")

    if args.scsplit is not None:
        scsplit_res = args.scsplit.split(':')
        scsplit_summary(scsplit_res, adata, mudata)
        print("scSplit result found")

    if args.souporcell is not None:
        souporcell_res = args.souporcell.split(':')
        souporcell_summary(souporcell_res, adata, mudata)
        print("Souporcell result found")

    # Read and combine assignment files
    assignment = [file for file in os.listdir("genetic_summary") if file.endswith("_assignment.csv")]
    assignment_all = pd.read_csv(os.path.join("genetic_summary", assignment[0]))

    if len(assignment) > 1:
        for df in assignment[1:]:
            df = pd.read_csv(os.path.join("genetic_summary", df))
            assignment_all = pd.merge(assignment_all, df, on='Barcode', how='outer')
    assignment_all.to_csv("genetic_summary/genetic_assignment_all.csv", index=False)

    # Read and combine classification files
    classification = [file for file in os.listdir("genetic_summary") if file.endswith("_classification.csv")]
    classification_all = pd.read_csv(os.path.join("genetic_summary", classification[0]))

    if len(classification) > 1:
        for df in classification[1:]:
            df = pd.read_csv(os.path.join("genetic_summary", df))
            classification_all = pd.merge(classification_all, df, on='Barcode', how='outer')
    classification_all.to_csv("genetic_summary/genetic_classification_all.csv", index=False)