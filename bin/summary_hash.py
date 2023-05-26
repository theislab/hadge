#!/usr/bin/env python
import pandas as pd
import os
import scanpy as sc
import argparse
import numpy as np
import muon as mu

parser = argparse.ArgumentParser(description="Parameters for summary process")
parser.add_argument("--demuxem", help="Folder containing output files of demuxem", default=None)
parser.add_argument("--htodemux", help="Folder containing output files of htodemux", default=None)
parser.add_argument("--multiseq", help="Folder containing output files of multiseq", default=None)
parser.add_argument("--hashsolo", help="Folder containing output files of hashsolo", default=None)
parser.add_argument("--solo", help="Folder containing output files of solo", default=None)
parser.add_argument("--hashedDrops", help="Folder containing output files of hashedDrops", default=None)
parser.add_argument("--demuxmix", help="Folder containing output files of Demuxmix", default=None)
parser.add_argument("--bff", help="Folder containing output files of BFF", default=None)
parser.add_argument("--gmm_demux", help="Folder containing output files of GMM-Demux", default=None)
parser.add_argument("--generate_anndata", help="Generate anndata", action='store_true')
parser.add_argument("--generate_mudata", help="Generate mudata", action='store_true')
parser.add_argument("--read_rna_mtx", help="10x-Genomics-formatted mtx directory for gene expression", default=None)
parser.add_argument("--read_hto_mtx", help="10x-Genomics-formatted mtx directory for HTO expression", default=None)

args = parser.parse_args()

def demuxem_summary(demuxem_res, raw_adata, raw_mudata):
    assign = []
    classi = []
    params = []
    for x in demuxem_res:
        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_obs.csv")][0])
        obs_res = pd.read_csv(obs_res_dir)
        obs_res.rename(columns={obs_res.columns[0]: "Barcode"}, inplace=True)
        demuxem_assign = obs_res[["Barcode", "assignment"]]
        demuxem_assign.rename(columns={"assignment": os.path.basename(x)}, inplace=True)
        demuxem_assign["Barcode"] = demuxem_assign["Barcode"].astype(str) + "-1"
        demuxem_assign.index = demuxem_assign.Barcode
        demuxem_assign = demuxem_assign.drop(columns=['Barcode'])
        assign.append(demuxem_assign)
        
        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(demuxem_assign, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")
        
        # TODO: add mudata

        demuxem_classi = obs_res[["Barcode", "demux_type"]]
        demuxem_classi.rename(columns={"demux_type": os.path.basename(x)}, inplace=True)
        demuxem_classi = demuxem_classi.replace("unknown", "negative")
        demuxem_classi["Barcode"] = demuxem_classi["Barcode"].astype(str) + "-1"
        demuxem_classi.index = demuxem_classi.Barcode
        demuxem_classi = demuxem_classi.drop(columns=['Barcode'])
        classi.append(demuxem_classi)

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
        #params_res.rename(columns={params_res.columns[1]: os.path.basename(x)}, inplace=True)
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
        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_res.csv")][0])
        obs_res = pd.read_csv(obs_res_dir, index_col=0)
        obs_res.index.name='Barcode'
        hashsolo_assign = obs_res[["Classification"]]
        hashsolo_assign.columns = [os.path.basename(x)]
        hashsolo_assign.replace({"Doublet": "doublet", "Negative": "negative"}, inplace=True)
        assign.append(hashsolo_assign)
        
        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(hashsolo_assign, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")
        
        hashsolo_classi = obs_res[["most_likely_hypothesis"]]
        hashsolo_classi["most_likely_hypothesis"] = hashsolo_classi["most_likely_hypothesis"].replace({0: "negative", 1: "singlet", 2: "doublet"})
        hashsolo_classi.columns = [os.path.basename(x)]
        classi.append(hashsolo_classi)
    
        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
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
        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_res.csv")][0])
        obs_res = pd.read_csv(obs_res_dir)
        
        obs_res["Classification"] = np.where(obs_res["Confident"], "singlet", np.where(obs_res["Doublet"], "doublet", "negative"))
        obs_res["Best"] = np.where(~obs_res["Classification"].isin(["doublet", "negative"]), obs_res["Best"], obs_res["Classification"])
        obs_res.rename(columns={obs_res.columns[0]: "Barcode"}, inplace=True)
        
        hasheddrops_res = obs_res[["Barcode", "Best"]]
        hasheddrops_res.rename(columns={"Best": os.path.basename(x)}, inplace=True)
        assign.append(hasheddrops_res)

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(hasheddrops_res, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")
        
        hasheddrops_classi = obs_res[["Barcode", "Classification"]]
        hasheddrops_classi.rename(columns={"Classification": os.path.basename(x)}, inplace=True)
        classi.append(hasheddrops_classi)
    
        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
        params_res = pd.read_csv(params_dir, usecols=[1, 2], keep_default_na=False, index_col=0)
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1).reset_index(drop=True)
    assign.to_csv("hash_summary/hasheddrops_assignment.csv", index=False, quoting=False)
    
    classi = pd.concat(classi, axis=1).reset_index(drop=True)
    classi.to_csv("hash_summary/hasheddrops_classification.csv", index=False, quoting=False)
    
    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/hasheddrops_params.csv")

def multiseq_summary(multiseq_res, raw_adata, raw_mudata):
    assign = []
    params = []
    for x in multiseq_res:
        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_res.csv")][0])
        multiseq_assign = pd.read_csv(obs_res_dir)
        multiseq_assign.columns = ["Barcode", os.path.basename(x)]
        multiseq_assign.index = multiseq_assign.Barcode
        multiseq_assign = multiseq_assign.drop(columns=['Barcode'])
        multiseq_assign.replace("Doublet", "doublet", inplace=True)
        multiseq_assign.replace("Negative", "negative", inplace=True)
        assign.append(multiseq_assign)
        
        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(multiseq_assign, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
        params_res = pd.read_csv(params_dir, usecols=[1, 2], keep_default_na=False, index_col=0)
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
        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_assignment_htodemux.csv")][0])
        htodemux_assign = pd.read_csv(obs_res_dir)
        htodemux_assign.columns = ["Barcode", os.path.basename(x)]
        htodemux_assign.replace("Doublet", "doublet", inplace=True)
        htodemux_assign.replace("Negative", "negative", inplace=True)
        htodemux_assign.index = htodemux_assign.Barcode
        htodemux_assign = htodemux_assign.drop(columns=['Barcode'])
        assign.append(htodemux_assign)

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(htodemux_assign, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")

        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_classification_htodemux.csv")][0])
        htodemux_classi = pd.read_csv(obs_res_dir)
        htodemux_classi.columns = ["Barcode", os.path.basename(x)]
        htodemux_classi.replace("Singlet", "singlet", inplace=True)
        htodemux_classi.replace("Doublet", "doublet", inplace=True)
        htodemux_classi.replace("Negative", "negative", inplace=True)
        htodemux_classi.index = htodemux_classi.Barcode
        htodemux_classi = htodemux_classi.drop(columns=['Barcode'])
        classi.append(htodemux_classi)


        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename == "params.csv"][0])
        params_res = pd.read_csv(params_dir, usecols=[1, 2], keep_default_na=False, index_col=0)     
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    assign = pd.concat(assign, axis=1)
    assign.to_csv("hash_summary/htodemux_assignment.csv", quoting=False)

    classi = pd.concat(classi, axis=1)
    classi.to_csv("hash_summary/htodemux_classification.csv", quoting=False)
    
    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/htodemux_params.csv")
        
def solo_summary(solo_res, raw_adata, raw_mudata):
    classi = []
    params = []
    
    for x in solo_res:
        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_res.csv")][0])
        solo_classi = pd.read_csv(obs_res_dir)
        solo_classi.columns = ["Barcode", os.path.basename(x)]
        solo_classi["Barcode"] = solo_classi["Barcode"].str.replace("-0", "")
        solo_classi.index = solo_classi.Barcode
        solo_classi = solo_classi.drop(columns=['Barcode'])
        classi.append(solo_classi)

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(solo_classi, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'droplet'}, inplace=True)
            adata.obs.droplet = adata.obs.droplet.fillna("negative")
            adata.obs.droplet = adata.obs.droplet.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
        params_res = pd.read_csv(params_dir, usecols=[0, 1], keep_default_na=False, index_col=0)
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)
    
    classi = pd.concat(classi, axis=1, join="outer")
    classi.to_csv("hash_summary/solo_classification.csv", quoting=False)
    
    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/solo_params.csv")


def demuxmix_summary(demuxmix_res,raw_adata, raw_mudata):
    classi = []
    assign = []

    for x in demuxmix_res:
        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_assignment_demuxmix.csv.csv")][0])
        
        demuxmix_asign = pd.read_csv(obs_res_dir)
        dt_assign = demuxmix_asign[["Barcode", "HTO"]]
        dt_assign.columns = ["Barcode", os.path.basename(x)]
        assign.append(dt_assign)

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(demuxmix_asign, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")

        demuxmix_classi = pd.read_csv(obs_res_dir)
        dt_classi = demuxmix_classi[["Barcode", "Classification"]]
        dt_classi.columns = ["Barcode", os.path.basename(x)]
        classi.append(dt_classi)

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename == "params.csv"][0])
        params_res = pd.read_csv(params_dir, usecols=[1, 2], keep_default_na=False, index_col=0)     
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    classi_df = pd.concat(classi, axis=1, join="outer")
    classi_df.to_csv("hash_summary/demuxmix_classification.csv", quoting=pd.QUOTE_NONE)
    
    assign_df = pd.concat(assign, axis=1, join="outer")
    assign_df.to_csv("hash_summary/demuxmix_assignment.csv", quoting=pd.QUOTE_NONE)

    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/demuxmix_params.csv")

def gmm_summary(gmmDemux_res,raw_adata, raw_mudata):
    classi = []
    assign = []
    params = []

    for x in gmmDemux_res:
        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("GMM_full.csv")][0])
        
        gmm_classi = pd.read_csv(obs_res_dir)
        
        #GMM full is the name given by GMM_demux per default to all results
        classification_config = os.path.join(x, "GMM_full.config")
        classif_file = pd.read_csv(classification_config)
        
        #Classification and Assigment come from the same file
        gmm_dt = pd.DataFrame(gmm_classi)
        classification_dt = pd.DataFrame(classif_file)
        
        classification_dt["V2"] = classification_dt["V2"].str.contains("-")
        classification_dt["Classification"] = classification_dt["V2"].apply(lambda x: "doublet" if x else "singlet")
        classification_dt.at[0, "Classification"] = "negative"
        gmm_dt["Classification"] = "negative"

        gmm_dt_assign = pd.DataFrame(gmm_classi)
        gmm_dt_assign["Assignment"] = result.map(classification_dt.set_index("V1")["V2"]).fillna(pd.NA)
        gmm_dt_assign = gmm_dt_assign[["V1", "Assignment"]]
        gmm_dt_assign.columns = ["Barcode", os.path.basename(x)]
        gmm_dt_assign["Barcode"] = gmm_dt_assign["Barcode"].str.replace("-1", "")
        assign.append(gmm_dt_assign)

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(gmm_dt_assign, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")

        
        result = gmm_dt["Cluster_id"].apply(lambda x: any(x == classification_dt["V1"]))
        gmm_dt["Classification"] = result.map(classification_dt.set_index("V1")["Classification"]).fillna(pd.NA)
        
        gmm_dt = gmm_dt[["V1", "Classification"]]
        gmm_dt.columns = ["Barcode", os.path.basename(x)]
        
        gmm_dt["Barcode"] = gmm_dt["Barcode"].str.replace("-1", "")
        
        classi.append(gmm_dt)

        params_dir = os.path.join(x, "params.csv")
        params_res = pd.read_csv(params_dir)
        params_res.columns = ["Argument", os.path.basename(x)]
        params.append(params_res)

    classi_df = pd.concat(classi, axis=1, join="outer")
    classi_df.to_csv("hash_summary/gmm_classification.csv", index=False, quoting=pd.QUOTE_NONE)
    
    
    assign_df = pd.concat(assign, axis=1, join="outer")
    assign_df.to_csv("hash_summary/GMM_assignment.csv", index=False, quoting=pd.QUOTE_NONE)
    
    
    params_df = pd.concat(params, axis=1, join="outer")
    params_df.to_csv("hash_summary/GMM_params.csv", index=False, quoting=pd.QUOTE_NONE)       

def bff_summary(bff_res,raw_adata, raw_mudata):
    classi = []
    assign = []

    for x in bff_res:
        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_bff.csv")][0])

        bff_assign = pd.read_csv(obs_res_dir)
        data_bff = pd.DataFrame(bff_assign)
        dt_assign = data_bff[["cellbarcode", "consensuscall"]]
        dt_assign = dt_assign.rename(columns={"cellbarcode": "Barcode", "consensuscall": os.path.basename(x)})
        dt_assign.loc[dt_assign["consensuscall"] == "Doublet", "consensuscall"] = "doublet"
        assign.append(dt_assign)

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(bff_assign, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")

        bff_classi = pd.read_csv(obs_res_dir)
        data_bff = pd.DataFrame(bff_classi)
        dt_classi = data_bff[["cellbarcode", "consensuscall.global"]]
        dt_classi = dt_classi.rename(columns={"cellbarcode": "Barcode", "consensuscall.global": os.path.basename(x)})
        dt_classi.loc[dt_classi["consensuscall.global"] == "Singlet", "consensuscall.global"] = "singlet"
        dt_classi.loc[dt_classi["consensuscall.global"] == "Doublet", "consensuscall.global"] = "doublet"
        dt_classi.loc[dt_classi["consensuscall.global"] == "Negative", "consensuscall.global"] = "negative"
        classi.append(dt_classi)

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename == "params.csv"][0])
        params_res = pd.read_csv(params_dir, usecols=[1, 2], keep_default_na=False, index_col=0)     
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    classi_df = pd.concat(classi, axis=1, join="outer")
    classi_df.to_csv("hash_summary/bff_classification.csv", index=False, quoting=pd.QUOTE_NONE)
    
    assign_df = pd.concat(assign, axis=1, join="outer")
    assign_df.to_csv("hash_summary/bff_assignment.csv", index=False, quoting=pd.QUOTE_NONE)
    
    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary/demuxmix_params.csv")

if __name__ == '__main__':
    adata = None
    mudata = None

    os.mkdir("hash_summary")
    os.mkdir("hash_summary/adata")
    os.mkdir("hash_summary/mudata")

    if args.generate_anndata is True:
        adata = sc.read_10x_mtx(args.read_rna_mtx)
    
    if args.generate_mudata is True:
        # TODO
        pass

    if args.hashedDrops is not None:
        hashedDrops_res = args.hashedDrops.split(':')
        hasheddrops_summary(hashedDrops_res, adata, mudata)
        print("hashedDrops result found")

    if args.demuxem is not None:
        demuxem_res = args.demuxem.split(':')
        demuxem_summary(demuxem_res, adata, mudata)
        print("DemuxEM result found")

    if args.hashsolo is not None:
        hashsolo_res = args.hashsolo.split(':')
        hashsolo_summary(hashsolo_res, adata, mudata)
        print("HashSolo result found")

    if args.multiseq is not None:
        multiseq_res = args.multiseq.split(':')
        multiseq_summary(multiseq_res, adata, mudata)
        print("MultiSeqDemux result found")

    if args.htodemux is not None:
        htodemux_res = args.htodemux.split(':')
        htodemux_summary(htodemux_res, adata, mudata)
        print("HTODemux result found")

    if args.solo is not None:
        solo_res = args.solo.split(':')
        solo_summary(solo_res, adata, mudata)
        print("solo result found")
    
    if args.demuxmix is not None:
        demuxmix_res = args.solo.split(':')
        demuxmix_summary(demuxmix_res, adata, mudata)
        print("Demuxmix result found")

    if args.gmm_demux is not None:
        gmmDemux_res = args.solo.split(':')
        gmm_summary(gmmDemux_res, adata, mudata)
        print("BFF result found")

    if args.bff is not None:
        bff_res = args.solo.split(':')
        bff_summary(bff_res, adata, mudata)
        print("GMM-Demux result found")


    # Read and combine assignment files
    assignment = [file for file in os.listdir("hash_summary") if file.endswith("_assignment.csv")]
    assignment_all = pd.read_csv(os.path.join("hash_summary", assignment[0]))

    if len(assignment) > 1:
        for df in assignment[1:]:
            df = pd.read_csv(os.path.join("hash_summary",df))
            assignment_all = pd.merge(assignment_all, df, on='Barcode', how='outer')
    assignment_all.to_csv("hash_summary/hashing_assignment_all.csv", index=False)

    # Read and combine classification files
    classification = [file for file in os.listdir("hash_summary") if file.endswith("_classification.csv")]
    classification_all = pd.read_csv(os.path.join("hash_summary", classification[0]))

    if len(classification) > 1:
        for df in classification[1:]:
            df = pd.read_csv(os.path.join("hash_summary",df))
            classification_all = pd.merge(classification_all, df, on='Barcode', how='outer')
    classification_all.to_csv("hash_summary/hashing_classification_all.csv", index=False)