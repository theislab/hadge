#!/usr/bin/env python
import pegasus as pg
import demuxEM
import numpy as np
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Parser for DemuxEM - Demultiplexing')
parser.add_argument('--rna_matrix_dir', help= 'cellranger output folder which contains raw RNA count matrix in mtx format.')
parser.add_argument('--hto_matrix_dir', help= 'cellranger output folder which contains raw HTO (antibody tag) count matrix in mtx format.')
parser.add_argument('--randomState', help='Random seed set for reproducing results.', type=int, default=0)
parser.add_argument('--min_signal', help='Any cell/nucleus with less than min_signal hashtags from the signal will be marked as unknown.', type=float, default=10.0)
parser.add_argument('--min_num_genes', help='We only demultiplex cells/nuclei with at least <number> of expressed genes.', type=int, default=100)
parser.add_argument('--min_num_umis', help='We only demultiplex cells/nuclei with at least <number> of UMIs.', type=int, default=100)
parser.add_argument('--alpha', help='The Dirichlet prior concentration parameter (alpha) on samples. An alpha value < 1.0 will make the prior sparse.', type=float, default=0.0)
parser.add_argument('--alpha_noise', help='The Dirichlet prior concenration parameter on the background noise.', type=float, default=1.0)
parser.add_argument('--tol', help='Threshold used for the EM convergence.', type=float, default=1e-6)
parser.add_argument('--n_threads', help='Number of threads to use. Must be a positive integer.', type=int, default=1)
parser.add_argument('--generateGenderPlot', help='Generate violin plots using gender-specific genes (e.g. Xist). <gene> is a comma-separated list of gene names.', default='')
parser.add_argument('--objectOutDemuxem', help='Output name of demultiplexing results. All outputs will use it as the prefix.', default="demuxem_res")
parser.add_argument('--outputdir', help='Output directory')

args = parser.parse_args()
param_list = [['rna_matrix_dir', args.rna_matrix_dir], ['hto_matrix_dir', args.hto_matrix_dir], ['randomState', args.randomState], ['min_signal', args.min_signal], ['min_num_genes', args.min_num_genes], ['min_num_umis', args.min_num_umis], ['alpha', args.alpha], ['alpha_noise', args.alpha_noise], ['tol', args.tol], ['n_threads', args.n_threads], ['generateGenderPlot', args.generateGenderPlot]]
 
param_df = pd.DataFrame(param_list, columns=['Argument', 'Value'])

if __name__ == '__main__':
    output_name = args.outputdir + "/" + args.objectOutDemuxem
    # load input rna data
    data = pg.read_input(args.rna_matrix_dir, modality="rna")
    data.subset_data(modality_subset=['rna'])
    data.concat_data() # in case of multi-organism mixing data
    # load input hashing data
    data.update(pg.read_input(args.hto_matrix_dir, modality="hashing"))
    # Extract rna and hashing data
    rna_data = data.get_data(modality="rna")
    hashing_data = data.get_data(modality="hashing")
    # Filter the RNA matrix
    rna_data.obs["n_genes"] = rna_data.X.getnnz(axis=1)
    rna_data.obs["n_counts"] = rna_data.X.sum(axis=1).A1
    obs_index = np.logical_and.reduce(
        (
            rna_data.obs["n_genes"] >= args.min_num_genes,
            rna_data.obs["n_counts"] >= args.min_num_umis,
        )
    )
    rna_data._inplace_subset_obs(obs_index)
    # run demuxEM
    demuxEM.estimate_background_probs(hashing_data, random_state=args.randomState)
    demuxEM.demultiplex(rna_data, hashing_data, min_signal=args.min_signal, alpha=args.alpha, alpha_noise=args.alpha_noise, tol=args.tol, n_threads=args.n_threads)
    # annotate raw matrix with demuxEM results
    demux_results = demuxEM.attach_demux_results(args.rna_matrix_dir, rna_data)
    # generate plots
    demuxEM.plot_hto_hist(hashing_data, "hto_type", output_name + ".ambient_hashtag.hist.pdf", alpha=1.0)
    demuxEM.plot_bar(hashing_data.uns["background_probs"], hashing_data.var_names, "Sample ID",
         "Background probability", output_name + ".background_probabilities.bar.pdf",)
    demuxEM.plot_hto_hist(hashing_data, "rna_type", output_name + ".real_content.hist.pdf", alpha=0.5)
    demuxEM.plot_rna_hist(rna_data, hashing_data, output_name + ".rna_demux.hist.pdf")
    
    if len(args.generateGenderPlot) > 0:
        rna_data.matrices["raw.X"] = rna_data.X.copy()
        rna_data.as_float()
        scale = 1e5 / rna_data.X.sum(axis=1).A1
        rna_data.X.data *= np.repeat(scale, np.diff(data.X.indptr))
        rna_data.X.data = np.log1p(rna_data.X.data)

        for gene_name in args.generateGenderPlot:
            plot_gene_violin(
                rna_data,
                gene_name,
                "{output_name}.{gene_name}.violin.pdf".format(
                    output_name=output_name, gene_name=gene_name
                ),
                title="{gene_name}: a gender-specific gene".format(gene_name=gene_name),
            )
    # output results
    pg.write_output(demux_results, output_name + "_demux.zarr.zip")
    pg.write_output(data, output_name + ".out.demuxEM.zarr.zip")
    print("\nSummary statistics:")
    print("total\t{}".format(rna_data.shape[0]))
    for name, value in rna_data.obs["demux_type"].value_counts().iteritems():
        print("{}\t{}".format(name, value))
    summary = rna_data.obs["demux_type"].value_counts().rename_axis('classification').reset_index(name='counts')
    total = ["total", rna_data.shape[0]]
    summary.loc[len(summary)] = total
    summary.to_csv(output_name + "_summary.csv", index=False)
    param_df.to_csv(args.outputdir + "/params.csv", index=False)
    
    hashtags = hashing_data.var.index.tolist()
    toreplace = [ht for ht in rna_data.obs['assignment'].unique() if ht not in hashtags]
    rna_data.obs.replace(toreplace,'doublet', inplace=True)
    rna_data.obs.to_csv(output_name + "_obs.csv")
