#!/usr/bin/env python
import pandas as pd
import os
import scanpy as sc
import argparse
import numpy as np
from mudata import MuData

parser = argparse.ArgumentParser(description="Parameters for summary process")
parser.add_argument("--demuxem", help="Folder containing output files of demuxem", default=None)
parser.add_argument("--htodemux", help="Folder containing output files of htodemux", default=None)
parser.add_argument("--multiseq", help="Folder containing output files of multiseq", default=None)
parser.add_argument("--hashsolo", help="Folder containing output files of hashsolo", default=None)
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
        obs_res.loc[:, "barcodekey"] = obs_res["barcodekey"].apply(lambda x: x + "-1")
        demuxem_assign = obs_res[["barcodekey", "assignment"]]
        demuxem_assign.columns = ["Barcode",os.path.basename(x)]
        demuxem_assign.set_index("Barcode", inplace=True)
        assign.append(demuxem_assign)
        
        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(demuxem_assign, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")
        
        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata['rna'].obs = mudata['rna'].obs.merge(demuxem_assign, left_index=True, right_index=True, how='left')
            mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
            mudata.update()
            mudata.write("hash_summary/mudata/mudata_with_mudata_"+ os.path.basename(x)+".h5mu") 

        demuxem_classi = obs_res[["barcodekey", "demux_type"]]
        demuxem_classi.columns = ["Barcode", os.path.basename(x)]
        demuxem_classi = demuxem_classi.replace("unknown", "negative")
        demuxem_classi.set_index("Barcode", inplace=True)
        classi.append(demuxem_classi)

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
        params_res = pd.read_csv(params_dir, keep_default_na=False, index_col=0)
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
        hashsolo_assign = hashsolo_assign.replace({"Doublet": "doublet", "Negative": "negative"})
        assign.append(hashsolo_assign)
        
        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(hashsolo_assign, left_index=True, right_index=True, how='left')
            adata.obs = adata.obs.rename(columns={adata.obs.columns[0]: 'donor'})
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")
        
        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata['rna'].obs = mudata['rna'].obs.merge(hashsolo_assign, left_index=True, right_index=True, how='left')
            mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
            mudata.update()
            mudata.write("hash_summary/mudata/mudata_with_mudata_"+ os.path.basename(x)+".h5mu") 

        hashsolo_classi = obs_res[["most_likely_hypothesis"]]
        hashsolo_classi.loc[:, "most_likely_hypothesis"] = hashsolo_classi["most_likely_hypothesis"].replace({0.0: "negative", 1.0: "singlet", 2.0: "doublet"})
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
        hasheddrops_res = hasheddrops_res.rename(columns={"Best": os.path.basename(x)})
        assign.append(hasheddrops_res)

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(hasheddrops_res, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")
        
        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata['rna'].obs = mudata['rna'].obs.merge(hasheddrops_res, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
            mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
            mudata.update()
            
            mudata.write("hash_summary/mudata/mudata_with_mudata_"+ os.path.basename(x)+".h5mu") 

        hasheddrops_classi = obs_res[["Barcode", "Classification"]]
        hasheddrops_classi = hasheddrops_classi.rename(columns={"Classification": os.path.basename(x)})
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
        multiseq_assign.set_index("Barcode", inplace=True)
        multiseq_assign.replace({"Doublet": "doublet", "Negative": "negative"}, inplace=True)
        assign.append(multiseq_assign)
        
        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(multiseq_assign, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata['rna'].obs = mudata['rna'].obs.merge(multiseq_assign, left_index=True, right_index=True, how='left')
            mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
            mudata.update()
            mudata.write("hash_summary/mudata/mudata_with_mudata_"+ os.path.basename(x)+".h5mu") 

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("params.csv")][0])
        params_res = pd.read_csv(params_dir, usecols=[1, 2], keep_default_na=False, index_col=0)
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)
    
    assign = pd.concat(assign, axis=1)
    assign.to_csv("hash_summary/multiseq_assignment.csv", quoting=False)

    classi = assign.copy()
    classi[~classi.isin(["doublet", "negative"])] = "singlet"
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
        htodemux_assign.replace({"Doublet": "doublet", "Negative": "negative"}, inplace=True)
        htodemux_assign.set_index("Barcode", inplace=True)
        assign.append(htodemux_assign)

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(htodemux_assign, left_index=True, right_index=True, how='left')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")

        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata['rna'].obs = mudata['rna'].obs.merge(htodemux_assign, left_index=True, right_index=True, how='left')
            mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
            mudata.update()
            mudata.write("hash_summary/mudata/mudata_with_mudata_"+ os.path.basename(x)+".h5mu") 

        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_classification_htodemux.csv")][0])
        htodemux_classi = pd.read_csv(obs_res_dir)
        htodemux_classi.columns = ["Barcode", os.path.basename(x)]
        htodemux_classi.replace({"Doublet": "doublet", "Negative": "negative", "Singlet": "singlet"}, inplace=True)
        htodemux_classi.set_index("Barcode", inplace=True)
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
        
def demuxmix_summary(demuxmix_res,raw_adata, raw_mudata):
    classi = []
    assign = []
    params = []
    files_in_folder = [item for item in demuxmix_res if os.path.isfile(os.path.join(demuxmix_res, item))]

    if len(files_in_folder) > 0:
        for x in demuxmix_res:
            obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_assignment_demuxmix.csv")][0])
            demuxmix_asign = pd.read_csv(obs_res_dir)
            if demuxmix_asign.empty:
                #no results create empty dataframe for empty col
                column_names = ['Barcode', os.path.basename(x)]
                # Create an empty dataframe with only column names
                df = pd.DataFrame(columns=column_names)
                classi.append(df)
                assign.append(df)
            else:
                dt_assign = demuxmix_asign[["Barcode", "HTO"]]
                dt_assign.columns = ["Barcode", os.path.basename(x)]
                assign.append(dt_assign)

                if raw_adata is not None:
                    adata = raw_adata.copy()
                    adata.obs = adata.obs.merge(demuxmix_asign, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
                    adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
                    adata.obs.donor = adata.obs.donor.fillna("negative")
                    adata.obs.donor = adata.obs.donor.astype(str)
                    adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")

                if raw_mudata is not None:
                    mudata = raw_mudata.copy()
                    mudata['rna'].obs = mudata['rna'].obs.merge(demuxmix_asign, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
                    mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
                    mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
                    mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
                    mudata.update()
                    mudata.write("hash_summary/mudata/mudata_with_mudata_"+ os.path.basename(x)+".h5mu") 

                demuxmix_classi = pd.read_csv(obs_res_dir)
                dt_classi = demuxmix_classi[["Barcode", "Classification"]]
                dt_classi.columns = ["Barcode", os.path.basename(x)]
                classi.append(dt_classi) 

                params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename == "params.csv"][0])
                params_res = pd.read_csv(params_dir, usecols=[1, 2], keep_default_na=False, index_col=0)     
                params_res.columns = [os.path.basename(x)]
                params.append(params_res)

        classi_df = pd.concat(classi, axis=1, join="outer")
        classi_df.to_csv("hash_summary" + "/demuxmix_classification.csv",index=False)
        
        assign_df = pd.concat(assign, axis=1, join="outer")
        assign_df.to_csv("hash_summary" + "/demuxmix_assignment.csv",index=False)

        params = pd.concat(params, axis=1)
        params.to_csv("hash_summary" + "/demuxmix_params.csv",index=False)
    else:
        print("No results found for Demuxmix")

def gmm_summary(gmmDemux_res,raw_adata, raw_mudata):
    classi = []
    assign = []
    params = []
    for x in gmmDemux_res:
        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("GMM_full.csv")][0])
        
        gmm_classi = pd.read_csv(obs_res_dir)
        
        #GMM full is the name given by GMM_demux per default to all results
        classification_config = os.path.join(x, "GMM_full.config")
        classif_file = pd.read_csv(classification_config,header=None)
        
        #Classification and Assigment come from the same file
        gmm_dt = pd.DataFrame(gmm_classi)
        classification_dt = pd.DataFrame(classif_file)
        #change column names
        classification_dt = classification_dt.rename(columns={0: "Cluster_id", 1: "assignment"})
        gmm_dt = gmm_dt.rename(columns={"Unnamed: 0": "Barcode"})
        #Create classification following the assignment found for the barcodes
        #we keep the original assigment and add a classification column
        classification_dt["assignment_binary"] = classification_dt["assignment"].str.contains("-")
        classification_dt["classification"] = classification_dt["assignment"].apply(lambda x: "doublet" if x else "singlet")
        classification_dt.at[0, "classification"] = "negative"
        #Compare classification guide file with classification found
        merged = pd.merge(classification_dt, gmm_dt, on='Cluster_id', how='left')

        gmm_dt['Classification'] = merged['classification']
        gmm_dt['Assignment'] = merged['assignment']
        
        #Assigment for GMM-Demux
        gmm_dt_assign = gmm_dt.drop(['Cluster_id','Confidence','Classification' ], axis=1)
        gmm_dt_assign.columns = ["Barcode", os.path.basename(x)]
        assign.append(gmm_dt_assign)

        if raw_adata is not None:
            adata = raw_adata.copy()
            adata.obs = adata.obs.merge(gmm_dt_assign, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
            adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
            adata.obs.donor = adata.obs.donor.fillna("negative")
            adata.obs.donor = adata.obs.donor.astype(str)
            adata.write("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")
        
        if raw_mudata is not None:
            mudata = raw_mudata.copy()
            mudata['rna'].obs = mudata['rna'].obs.merge(gmm_dt_assign, left_index=True, right_on='Barcode', how='left').set_index('Barcode')
            mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
            mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
            mudata.update()
            mudata.write("hash_summary/mudata/mudata_with_mudata_"+ os.path.basename(x)+".h5mu") 

        #Classification for GMM-Demux
        gmm_dt_classi = gmm_dt.drop(['Cluster_id','Confidence','Assignment' ], axis=1)
        gmm_dt_classi.columns =["Barcode", os.path.basename(x)]
        classi.append(gmm_dt_classi)

        params_dir = os.path.join(x, "params.csv")
        params_res = pd.read_csv(params_dir,index_col=False)
        params_res.columns = ["Argument", os.path.basename(x)]
        params.append(params_res)

    classi_df = pd.concat(classi, axis=1, join="outer")
    classi_df.to_csv("hash_summary"  +"/GMM_classification.csv", index=False)
    
    assign_df = pd.concat(assign, axis=1, join="outer")
    assign_df.to_csv("hash_summary"  +"/GMM_assignment.csv", index=False, sep=",")
    
    params_df = pd.concat(params, axis=1, join="outer")
    params_df.to_csv("hash_summary"  +"/GMM_params.csv", index=False)    

def bff_summary(bff_res,raw_adata, raw_mudata):
    classi = []
    assign = []
    params = []

    for x in bff_res:
        obs_res_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename.endswith("_bff.csv")][0])

        bff_assign = pd.read_csv(obs_res_dir)
        data_bff = pd.DataFrame(bff_assign)
        if data_bff.empty:
            #no results create empty dataframe for empty col
            column_names = ['Barcode', os.path.basename(x)]
            # Create an empty dataframe with only column names
            df = pd.DataFrame(columns=column_names)
            classi.append(df)
            assign.append(df)
        else:
            #df contain data and we save it in the same way
            dt_assign = data_bff.copy()
            #check if the columns contain results from both bff_s or only one
            column_names = ["Unnamed: 0","bff_raw","bff_cluster","consensuscall.global"]
            for column in column_names:
                if column in dt_assign.columns:
                    dt_assign = dt_assign.drop([column], axis=1)
            dt_assign.loc[dt_assign["consensuscall"] == "Singlet", "consensuscall"] = "singlet"
            dt_assign.loc[dt_assign["consensuscall"] == "Doublet", "consensuscall"] = "doublet"
            dt_assign = dt_assign.rename(columns={"cellbarcode": "Barcode", "consensuscall": os.path.basename(x)})
            assign.append(dt_assign)
            if raw_adata is not None:
                adata = raw_adata.copy()

                adata.obs = adata.obs.merge(dt_assign, left_index=True, right_index=True, how='left')

                adata.obs.rename(columns={adata.obs.columns[0]: 'donor'}, inplace=True)
                adata.obs.donor = adata.obs.donor.fillna("negative")
                adata.obs.donor = adata.obs.donor.astype(str)
                adata.write_h5ad("hash_summary/adata/adata_with_"+os.path.basename(x)+".h5ad")
                

            if raw_mudata is not None:
                mudata = raw_mudata.copy()

                mudata['rna'].obs = mudata['rna'].obs.merge(dt_assign, left_index=True, right_index=True, how='left')
                mudata['rna'].obs.rename(columns={mudata['rna'].obs.columns[0]: 'donor'}, inplace=True)
                mudata['rna'].obs.donor = mudata['rna'].obs.donor.fillna("negative")
                mudata['rna'].obs.donor = mudata['rna'].obs.donor.astype(str)
                mudata.update()
                mudata.write("hash_summary/mudata/mudata_with_mudata_"+ os.path.basename(x)+".h5mu") 


            dt_classi = data_bff.copy()
            column_names_class = ["bff_raw","bff_cluster","consensuscall"]
            for column in column_names_class:
                if column in dt_assign.columns:
                    dt_classi = dt_classi.drop([column], axis=1)
            dt_classi.loc[dt_classi["consensuscall.global"] == "Singlet", "consensuscall.global"] = "singlet"
            dt_classi.loc[dt_classi["consensuscall.global"] == "Doublet", "consensuscall.global"] = "doublet"
            dt_classi.loc[dt_classi["consensuscall.global"] == "Negative", "consensuscall.global"] = "negative"
            dt_classi = dt_classi.rename(columns={"cellbarcode": "Barcode", "consensuscall.global": os.path.basename(x)})
            classi.append(dt_classi)

        params_dir = os.path.join(x, [filename for filename in os.listdir(x) if filename == "params.csv"][0])
        params_res = pd.read_csv(params_dir, usecols=[1, 2], keep_default_na=False, index_col=0)     
        params_res.columns = [os.path.basename(x)]
        params.append(params_res)

    classi_df = pd.concat(classi, axis=1, join="outer")
    classi_df.to_csv("hash_summary" +"/bff_classification.csv", index=False)
        
    assign_df = pd.concat(assign, axis=1, join="outer")
    assign_df.to_csv("hash_summary"  +"/bff_assignment.csv", index=False)
        
    params = pd.concat(params, axis=1)
    params.to_csv("hash_summary"  +"/bff_params.csv")

if __name__ == '__main__':
    adata = None
    mudata = None

    if not os.path.exists("hash_summary"):
        os.mkdir("hash_summary")

    if args.generate_anndata is True:
        if not os.path.exists("hash_summary/adata"):
            os.mkdir("hash_summary/adata")
        adata = sc.read_10x_mtx(args.read_rna_mtx)

    if args.generate_mudata is True:
        if not os.path.exists("hash_summary/mudata"):
            os.mkdir("hash_summary/mudata")
        rna_data = sc.read_10x_mtx(args.read_rna_mtx)
        hto_data = sc.read_10x_mtx(args.read_hto_mtx, gex_only=False)
        mudata = MuData({"rna": rna_data, "hto": hto_data })
        
    if args.hashedDrops is not None:
        hashedDrops_res = args.hashedDrops.split(':')
        hasheddrops_summary(hashedDrops_res, adata, mudata)

    if args.demuxem is not None:
        demuxem_res = args.demuxem.split(':')
        demuxem_summary(demuxem_res, adata, mudata)

    if args.hashsolo is not None:
        hashsolo_res = args.hashsolo.split(':')
        hashsolo_summary(hashsolo_res, adata, mudata)

    if args.multiseq is not None:
        multiseq_res = args.multiseq.split(':')
        multiseq_summary(multiseq_res, adata, mudata)

    if args.htodemux is not None:
        htodemux_res = args.htodemux.split(':')
        htodemux_summary(htodemux_res, adata, mudata)
    
    if args.demuxmix is not None:
        demuxmix_res = args.demuxmix.split(':')
        demuxmix_summary(demuxmix_res, adata, mudata)

    if args.gmm_demux is not None:
        gmmDemux_res = args.gmm_demux.split(':')
        gmm_summary(gmmDemux_res, adata, mudata)

    if args.bff is not None:
        bff_res = args.bff.split(':')
        bff_summary(bff_res, adata, mudata)

    # Read and combine assignment files
    assignment = [file for file in os.listdir("hash_summary") if file.endswith("_assignment.csv")]
    assignment_all = pd.read_csv(os.path.join("hash_summary", assignment[0]))
   
    if len(assignment) > 1:
        for df in assignment[1:]:
            df = pd.read_csv(os.path.join("hash_summary", df))
            assignment_all = pd.merge(assignment_all, df, on='Barcode', how='outer')
    assignment_all.to_csv("hash_summary/hashing_assignment_all.csv", index=False)

    # Read and combine classification files
    classification = [file for file in os.listdir("hash_summary") if file.endswith("_classification.csv")]
    classification_all = pd.read_csv(os.path.join("hash_summary", classification[0]))

    if len(classification) > 1:
        for df in classification[1:]:
            df = pd.read_csv(os.path.join("hash_summary", df))
            classification_all = pd.merge(classification_all, df, on='Barcode', how='outer')
    classification_all.to_csv("hash_summary/hashing_classification_all.csv", index=False)