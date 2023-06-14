# Usage

## **Input data preparation**

The input data depends heavily on the deconvolution tools. In the following table, you will find the minimal input data required by different tools.

### Genotype-based deconvolution methods:

| Deconvolution methods | Input data                                                                           |
| --------------------- | ------------------------------------------------------------------------------------ |
| Demuxlet              | - Alignment (BAM)<br>- Barcode (TSV)<br>- Genotype reference (VCF)                   |
| Freemuxlet            | - Alignment (BAM)<br>- Barcode (TSV)<br>- Genotypes from referenced population (VCF) |
| Vireo                 | - Cell genotype (VCF or cellSNP folder)                                              |
| Souporcell            | - Alignment (BAM)<br>- Barcode (TSV)<br>- Reference genome (FASTA)                   |
| scSplit               | - Alignment (BAM)<br>- Barcode (TSV)<br>- Genotypes from referenced population (VCF) |
|                       |                                                                                      |

You may see that some tools share some input data in common, so we set only one parameter for the same input for benchmarking.

| Input data                                 | Parameter                              |
| ------------------------------------------ | -------------------------------------- |
| Alignment (BAM)                            | `params.bam`<br>`params.bai`           |
| Barcode (TSV)                              | `params.barcodes`                      |
| Genotype reference (VCF)                   | `params.vcf_donor`                     |
| Genotypes from referenced population (VCF) | `params.vcf_mixed`                     |
| Reference genome (FASTA)                   | `params.fasta`<br>`params.fasta_index` |
| Cell genotype (VCF or cellSNP folder)      | `params.celldata`                      |

#### Pre-processing

In case you want to perform genotype-based deconvolution on pre-processed data, we provide a process in concordance with [the instruction of scSplit](https://github.com/jon-xu/scSplit). It only requires the Alignment (BAM) file as input. To specify which method is performed on the pre-processed data : set `[method]_preprocess = True`.

#### Variant calling

In case you don't have any cell genotypes or variants called from mixed samples yet, we provide two processes for variant calling.
| Variant calling methods | Input data | Parameter | Output |
|------------------------- |------------------------------------------------------------- |-------------------------------------------------------------------------- |----------------------------- |
| freebayes | - Alignment (BAM)<br>- Reference genome (FASTA) | `params.bam`<br>`params.bai`<br>`params.fasta`<br>`params.fasta_index` | Variants from mixed samples |
| cellsnp-lite | - Alignment (BAM)<br>- Barcode (TSV)<br>- Common SNPs (VCF) | `params.bam`<br>`params.bai`<br>`params.barcodes`<br>`params.regionsVCF` | Cell genotypes |

You can have following options for `scsplit_variant`:

- `True`: activate freebayes
- Otherwise: inactivate variant calling, get the input data from `params.vcf_mixed`

You can have following options for `vireo_variant`:

- `True`: activate cellsnp
- Otherwise: inactivate variant calling, get the input data from `params.celldata`

#### Common variants

When running genotype-based deconvolution methods without genotype reference, you may need common variants from the popultion. Here we collect different sources of common variants for GRCh38 recommended by different methods.

| Method       | Paramter                   | Source                                                                    |
| ------------ | -------------------------- | ------------------------------------------------------------------------- |
| scSplit      | common_variants_scSplit    | https://melbourne.figshare.com/articles/dataset/Common_SNVS_hg38/17032163 |
| Souporcell   | common_variants_souporcell | https://github.com/wheaton5/souporcell                                    |
| Freemuxlet   | common_variants_freemuxlet | https://sourceforge.net/projects/cellsnp/files/SNPlist/                   |
| cellSNP-lite | common_variants_cellsnp    | https://sourceforge.net/projects/cellsnp/files/SNPlist/                   |

### Hashing-based deconvolution workflow

| Deconvolution method | Input data                                                                                                         | Parameter                                                      |
| -------------------- | ------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------- |
| HTODemux             | - Seurat object with both UMI and hashing count matrix (RDS)                                                       | `params.rna_matrix_htodemux` <br> `params.hto_matrix_htodemux` |
| Multiseq             | - Seurat object with both UMI and hashing count matrix (RDS)                                                       | `params.rna_matrix_multiseq` <br> `params.hto_matrix_multiseq` |
| HashSolo             | - 10x mtx directory with hashing count matrix (H5)                                                                 | `params.hto_matrix_hashsolo` <br> `params.rna_matrix_hashsolo` |
| HashedDrops          | - 10x mtx directory with hashing count matrix (Directory)                                                          | `params.hto_matrix_hashedDrops`                                |
| Demuxem              | - 10x mtx directory with UMI count matrix (Directory)<br>- 10x mtx directory with hashing count matrix (Directory) | `params.hto_matrix_demuxem`<br>`params.rna_matrix_demuxem`     |

The parameters `params.[rna/hto]_matrix_[method]` is used to specify whether to use raw or filtered counts for each method. Similary as genotype-based deconvlution methods, hashing methods also utilize common input parameters to store count matrices for better control.

| Input data                     | Parameter                    |
| ------------------------------ | ---------------------------- |
| Raw scRNAseq count matrix      | `params.rna_matrix_raw`      |
| Filtered scRNAseq count matrix | `params.rna_matrix_filtered` |
| Raw HTO count matrix           | `params.hto_matrix_raw`      |
| Filtered HTO count matrix      | `params.hto_matrix_filtered` |

#### Pre-processing

Similar as in the genetic demultiplexing workflow, we provide a pre-processing step required before running HTODemux and Multiseq to load count matrices into a Seurat object. The input will be automatically loaded from the parameters mentioned above.

### **Running on multiple samples**

The pipeline is able to run on multiple samples. In this scenario, the shared parameters for input data are retrieved from a sample sheet using `params.multi_sample`, which is set to None by default. Along with the input data, the sample sheet should contain an additional column for unique sample IDs assigned to each sample. The remaining parameters for each process are specified in the nextflow.config file, just like when demultiplexing a single sample. However, there is a distinction between running on a single sample and running on multiple samples. When processing multiple samples, the pipeline only permits a single value for each process parameter, whereas in the case of a single sample, multiple values separated by commas are allowed. The sample sheet should have e.g. following columns depending on the methods you want to run:

- sampleId
- na_matrix_raw
- rna_matrix_filtered
- hto_matrix_raw
- hto_matrix_filtered
- bam
- bam_index
- barcodes
- nsample
- celldata
- vcf_mixed
- vcf_donor

### **scverse compatibility**

To ensure scverse compatibility, the pipeline provides the option to generate anndata or mudata specifeid by `params.generate_anndata`. If set to True, the pipeline will generate an AnnData object in the folder `[workflow]_summary/adata` during the summary process of two workflows. This object contains the scRNA-seq counts from `params.rna_matrix_filered` and stores the assignment of each demultiplexing method in the `assignment` column of `obs`. Additionlly, if `match_donor` is True, the pipeline also produces an AnnData object in the `data_output` folder which contains the assignment of the best-matched method pair after donor matching.

## **Pipeline configuration**

### **Conda environments:**

We provide a `environment.yml` file for each process. But you can also use local Conda environments to run a process:

```
// dont forget to enable conda
conda.enable = true
process {
    // Use Conda environment files
    withName:scSplit {
        conda = './conda/scsplit.yml'
    }
    // Use Conda package names
    withName:cellSNP {
        conda = 'bioconda::cellsnp-lite'
    }
    // Use existing Conda environments
    withName:summary {
        conda = '/path/to/an/existing/env/directory'
    }
}

```

### Containers:

Nextflow also supports a variety of container runtimes, e.g. Docker. To specify a different Docker image for each process:

```
process {
    withName:foo {
        container = 'image_name_1'
    }
    withName:bar {
        container = 'image_name_2'
    }
}
// do not forget to enable docker

docker.enabled = true

```

### Executor and resource specifications:

- The pipeline can be run either locally or on an HPC. You can set the executor by running the pipeline with `-profile standard` or `-profile cluster`. Of course, you can add other profiles if you want.
- Feel free to add other configurations, e.g. the number of CPUS, the memory allocation, etc. If you are new to Nextflow framework, please visit the [Nextlfow page](https://www.nextflow.io/docs/latest/config.html#).
- As default, the pipeline is run locally with the standard profile, where all processes annotated with the big_mem label are assigned 4 cpus and 16 Gb of memory.

```
profiles{
    standard {
        process {
            executor = 'local'
            withLabel: big_mem {
                cpus = 4
                memory = 16.GB
            }
            withLabel: small_mem {
                cpus = 2
                memory = 8.GB
            }
        }

    }

    cluster {
        process {
            executor = 'slurm'
             // queue = ...
            withLabel: big_mem {
                cpus = 32
                memory = 64.GB
            }
            withLabel: small_mem {
                cpus = 16
                memory = 32.GB
            }
        }
    }
}

```

## Parameters

### General

|                  |                                                                 |
| :--------------: | :-------------------------------------------------------------: |
|      outdir      |                Output directory of the pipeline                 |
|       mode       |  Mode of the pipeline: genetic, hashing, rescue or donor_match  |
| generate_anndata | Whether to generate anndata after demultiplexing. Default: True |
| generate_mudata  | Whether to generate mudata after demultiplexing. Default: False |

### Hashing-based: Preprocessing

|               |                                                                                                 |
| ------------- | ----------------------------------------------------------------------------------------------- |
| ndelim        | For the initial identity calss for each cell, delimiter for the cell's column name. Default: \_ |
| sel_method    | The selection method used to choose top variable features. Default: mean.var.plot               |
| n_features    | Number of features to be used when finding variable features. Default: 2000                     |
| assay         | Assay name for HTO modality. Default: HTO                                                       |
| norm_method   | Method for normalization of HTO data. Default: CLR                                              |
| margin        | If performing CLR normalization, normalize across features (1) or cells (2). Default: 2         |
| preprocessOut | Name of the output Seurat object. Default: preprocessed                                         |

### Hashing-based: HTODemux

|                     |                                                                                                                                |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| htodemux            | Whether to perform Multiseq. Default: True                                                                                     |
| rna_matrix_htodemux | Whether to use raw or filtered scRNA-seq count matrix. Default: filtered                                                       |
| hto_matrix_htodemux | Whether to use raw or filtered HTO count matrix. Default: filtered                                                             |
| assay               | Name of the hashtag assay. Default: HTO                                                                                        |
| quantile_htodemux   | The quantile of inferred 'negative' distribution for each hashtag, over which the cell is considered 'positive'. Default: 0.99 |
| kfunc               | Clustering function for initial hashtag grouping. Default: clara.                                                              |
| nstarts             | nstarts value for k-means clustering when kfunc=kmeans. Default: 100                                                           |
| nsamples            | Number of samples to be drawn from the dataset used for clustering when kfunc= clara. Default: 100                             |
| seed                | Sets the random seed. Default: 42                                                                                              |
| init                | Initial number of clusters for hashtags. Default: NULL, which means the # of hashtag oligo names + 1 to account for negatives. |
| objectOutHTO        | Name of the output Seurat object. Default: htodemux                                                                            |
| assignmentOutHTO    | Prefix of the output CSV files. Default: htodemux                                                                              |
| ridgePlot           | Whether to generate a ridge plot to visualize enrichment for all HTOs. Default: TRUE                                           |
| ridgeNCol           | Number of columns in the ridge plot. Default: 3                                                                                |
| featureScatter      | Whether to generate a scatter plot to visualize pairs of HTO signals. Default: FALSE                                           |
| scatterFeat1        | First feature to plot. Default: None                                                                                           |
| scatterFeat2        | Second feature to plot. Default: None                                                                                          |
| vlnplot             | Whether to generate a violin plot, e.g. to compare number of UMIs for singlets, doublets and negative cells. Default: TRUE     |
| vlnFeatures         | Features to plot. Default: nCount_RNA                                                                                          |
| vlnLog              | Whether to plot the feature axis on log scale. Default: TRUE                                                                   |
| tsne                | Whether to generate a 2D tSNE embedding for HTOs. Default: TRUE                                                                |
| tsneIdents          | Subset Seurat object based on identity class. Default: Negative                                                                |
| tsneInvert          | Whether to keep or remove the identity class. Default: TRUE                                                                    |
| tsneVerbose         | Whether to print the top genes associated with high/low loadings for the PCs when running PCA. Default: FALSE                  |
| tsneApprox          | Whether to use truncated singular value decomposition to approximate PCA. Default: FALSE                                       |
| tsneDimMax          | Number of dimensions to use as input features when running t-SNE dimensionality reduction. Default: 2                          |
| tsnePerplexity      | Perplexity when running t-SNE dimensionality reduction. Default: 100                                                           |
| heatmap             | Whether to generate an HTO heatmap. Default: TRUE                                                                              |
| heatmapNcells       | Number of cells to plot. Default: 5000                                                                                         |

### Hashing-based: Multiseq

|                     |                                                                                                         |
| ------------------- | ------------------------------------------------------------------------------------------------------- |
| multiseq            | Whether to perform Multiseq. Default: True                                                              |
| rna_matrix_multiseq | Whether to use raw or filtered scRNA-seq count matrix. Default: filtered                                |
| hto_matrix_multiseq | Whether to use raw or filtered HTO count matrix. Default: filtered                                      |
| assay               | Name of the hashtag assay, same as used for HTODemux. Default: HTO                                      |
| quantile_multi      | The quantile to use for classification. Default: 0.7                                                    |
| autoThresh          | Whether to perform automated threshold finding to define the best quantile. Default: TRUE               |
| maxiter             | nstarts value for k-means clustering when kfunc=kmeans. Default: 100                                    |
| qrangeFrom          | The minimal possible quantile value to try if autoThresh=TRUE. Default: 0.1                             |
| qrangeTo            | The minimal possible quantile value to try if autoThresh=TRUE. Default: 0.9                             |
| qrangeBy            | The constant difference of a range of possible quantile values to try if autoThresh=TRUE. Default: 0.05 |
| verbose_multiseq    | Wether to print the output. Default: TRUE                                                               |
| assignmentOutMulti  | Prefix of the output CSV files. Default: multiseq                                                       |
| objectOutMulti      | Name of the output Seurat object. Default: multiseq                                                     |

### Hashing-based: Solo

|                            |                                                                                                  |
| -------------------------- | ------------------------------------------------------------------------------------------------ |
| solo                       | Whether to perform Solo. Default: True                                                           |
| rna_matrix_solo            | Input folder to RNA expression matrix in 10x format.                                             |
| max_epochs                 | Number of epochs to train for. Default: 400                                                      |
| lr                         | Learning rate for optimization. Default: 0.001                                                   |
| train_size                 | Size of training set in the range between 0 and 1. Default: 0.9                                  |
| validation_size            | Size of the test set. Default: 0.1                                                               |
| batch_size                 | Minibatch size to use during training. Default: 128                                              |
| early_stopping             | Adds callback for early stopping on validation_loss. Default: True                               |
| early_stopping_patience    | Number of times early stopping metric can not improve over early_stopping_min_delta. Default: 30 |
| early_stopping_min_delta   | Threshold for counting an epoch towards patience train(). Default: 10                            |
| soft                       | Return probabilities instead of class label. Default: False                                      |
| include_simulated_doublets | Return probabilities for simulated doublets as well.                                             |
| assignmentOutSolo          | Prefix of the output CSV files. Default: solo_predict                                            |

### Hashing-based: HashSolo

|                          |                                                                                              |
| ------------------------ | -------------------------------------------------------------------------------------------- |
| hashsolo                 | Whether to perform HashSolo. Default: True                                                   |
| rna_matrix_hashsolo      | Whether to use raw or filtered scRNA-seq count matrix. Default: raw                          |
| hto_matrix_hashsolo      | Whether to use raw or filtered HTO count matrix if use_rna_data is set to True. Default: raw |
| priors_negative          | Prior for the negative hypothesis. Default: 1/3                                              |
| priors_singlet           | Prior for the singlet hypothesis. Default: 1/3                                               |
| priors_doublet           | Prior for the doublet hypothesis. Default: 1/3                                               |
| pre_existing_clusters    | Column in the input data for how to break up demultiplexing. Default: None                   |
| use_rna_data             | Whether to use RNA counts for deconvolution. Default: False                                  |
| number_of_noise_barcodes | Number of barcodes to use to create noise distribution. Default: None                        |
| assignmentOutHashSolo    | Prefix of the output CSV files. Default: hashsolo                                            |
| plotOutHashSolo          | Prefix of the output figures. Default: hashsolo                                              |

### Hashing-based: DemuxEm

|                      |                                                                                                                               |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| demuxem              | Whether to perform Demuxem. Default: True                                                                                     |
| rna_matrix_demuxem   | Whether to use raw or filtered scRNA-seq count matrix. Default: raw                                                           |
| hto_matrix_demuxem   | Whether to use raw or filtered HTO count matrix. Default: raw                                                                 |
| threads_demuxem      | Number of threads to use. Must be a positive integer. Default: 1                                                              |
| alpha_demuxem        | The Dirichlet prior concentration parameter (alpha) on samples. An alpha value < 1.0 will make the prior sparse. Default: 0.0 |
| alpha_noise          | The Dirichlet prior concenration parameter on the background noise. Default: 1.0                                              |
| min_num_genes        | Filter cells/nuclei with at least specified number of expressed genes. Default: 100                                           |
| min_num_umis         | Filter cells/nuclei with at least specified number of UMIs. Default: 100                                                      |
| min_signal           | Any cell/nucleus with less than min_signal hashtags from the signal will be marked as unknown. Default: 10                    |
| tol                  | Threshold used for the EM convergence. Default: 1e-6                                                                          |
| generate_gender_plot | Generate violin plots using gender-specific genes (e.g. Xist). Value is a comma-separated list of gene names. Default: None   |
| random_state         | Random seed set for reproducing results. Default: 0                                                                           |
| objectOutDemuxem     | Prefix of the output files. Default: demuxem_res                                                                              |

### Hashing-based: HashedDrops

|                          |                                                                                                                                                                                                            |
| ------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| hashedDrops              | Whether to perform hashedDrops. Default: True                                                                                                                                                              |
| hto_matrix_hashedDrops   | Whether to use raw or filtered HTO count matrix. Default: raw                                                                                                                                              |
| lower                    | The lower bound on the total UMI count, at or below which all barcodes are assumed to correspond to empty droplets. Default: 100                                                                           |
| niters                   | The number of iterations to use for the Monte Carlo p-value calculations. Default: 10000                                                                                                                   |
| testAmbient              | Whether results should be returned for barcodes with totals less than or equal to lower. Default: TRUE                                                                                                     |
| ignore_hashedDrops       | The lower bound on the total UMI count, at or below which barcodes will be ignored. Default: NULL                                                                                                          |
| alpha_hashedDrops        | The scaling parameter for the Dirichlet-multinomial sampling scheme. Default: NULL                                                                                                                         |
| round                    | Whether to check for non-integer values in m and, if present, round them for ambient profile estimation. Default: TRUE                                                                                     |
| byRank                   | If set, this is used to redefine lower and any specified value for lower is ignored. Default: NULL                                                                                                         |
| isCellFDR                | FDR Threshold to filter the cells for empty droplet detection. Default: 0.01                                                                                                                               |
| objectOutEmptyDrops      | Prefix of the emptyDroplets output RDS object. Default: emptyDroplets                                                                                                                                      |
| assignmentOutEmptyDrops  | Prefix of the emptyDroplets output CSV file. Default: emptyDroplets                                                                                                                                        |
| ambient                  | Whether to use the relative abundance of each HTO in the ambient solution from emptyDrops, set TRUE only when testAmbient=TRUE. Default: FALSE                                                             |
| minProp                  | The ambient profile when ambient=NULL. Default: 0.05                                                                                                                                                       |
| pseudoCount              | The minimum pseudo-count when computing logfold changes. Default: 5                                                                                                                                        |
| constantAmbient          | Whether a constant level of ambient contamination should be used to estimate LogFC2 for all cells. Default: FALSE                                                                                          |
| doubletNmads             | The number of median absolute deviations (MADs) to use to identify doublets. Default: 3                                                                                                                    |
| doubletMin               | The minimum threshold on the log-fold change to use to identify doublets. Default: 2                                                                                                                       |
| doubletMixture           | Wwhether to use a 2-component mixture model to identify doublets. Default: FALSE                                                                                                                           |
| confidentNmads           | The number of MADs to use to identify confidently assigned singlets. Default: 3                                                                                                                            |
| confidenMin              | The minimum threshold on the log-fold change to use to identify singlets. Default: 2                                                                                                                       |
| combinations             | An integer matrix specifying valid combinations of HTOs. Each row corresponds to a single sample and specifies the indices of rows in x corresponding to the HTOs used to label that sample. Default: NULL |
| objectOutHashedDrops     | Prefix of the hashedDrops output RDS object. Default: hashedDrops                                                                                                                                          |
| assignmentOutHashedDrops | Prefix of the hashedDrops output CSV file. Default: hashedDrops                                                                                                                                            |

### Genotype-based: Demuxlet and dsc-pileup

|                     |                                                                                                                                                              |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| demuxlet            | Whether to run Demuxlet. Default: False                                                                                                                      |
| demuxlet_preprocess | Whether to perform pre-processing on the input params.bam for demuxlet. True: Perform pre-processing. Otherwise pre-processing is not called. Default: False |
| bam                 | Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed.                                                                                          |
| bai                 | Index of Input SAM/BAM/CRAM file.                                                                                                                            |
| barcodes            | List of cell barcodes to consider.                                                                                                                           |
| tag_group           | Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB Default: CB                |
| tag_UMI             | Tag representing UMIs. For 10x genomiucs, use UB. Default: UB                                                                                                |
| sm                  | List of sample IDs to compare to. Default: None (use all)                                                                                                    |
| vcf_donor           | Input VCF/BCF file, containing GT, GP or PL for donors. It also requires the AC and AN field if plp_freemuxlet=True.                                         |
| sm_list             | File containing the list of sample IDs to compare. Default: None                                                                                             |
| sam_verbose         | Verbose message frequency for SAM/BAM/CRAM. Default: 1000000                                                                                                 |
| vcf_verbose         | Verbose message frequency for VCF/BCF. Default: 10000                                                                                                        |
| skip_umi            | Do not generate [prefix].umi.gz file, which stores the regions covered by each barcode/UMI pair. Default: False                                              |
| cap_BQ              | Maximum base quality (higher BQ will be capped). Default: 40                                                                                                 |
| min_BQ              | Minimum base quality to consider (lower BQ will be skipped). Default: 13                                                                                     |
| min_MQ              | Minimum mapping quality to consider (lower MQ will be ignored). Default: 20                                                                                  |
| min_TD              | Minimum distance to the tail (lower will be ignored). Default: 0                                                                                             |
| excl_flag           | SAM/BAM FLAGs to be excluded. Default: 3844                                                                                                                  |
| min_total           | Minimum number of total reads for a droplet/cell to be considered. Default: 0                                                                                |
| min_uniq            | Minimum number of unique reads (determined by UMI/SNP pair) for a droplet/cell to be considered. Default: 0                                                  |
| min_snp             | Minimum number of SNPs with coverage for a droplet/cell to be considered. Default: 0                                                                         |
| min_umi             | Minimum number of UMIs for a droplet/cell to be considered. Default: 0                                                                                       |
| plp                 | Whether to call dsc-pileup. If set True, dsc-pileup will be called. It set False, will use SAM file to call Demuxlet. Default: False                         |
| field               | FORMAT field to extract the genotype, likelihood, or posterior from. Default: GT                                                                             |
| geno_error_offset   | Offset of genotype error rate. [error] = [offset] + [1-offset]_[coeff]_[1-r2]. Default: 0.1                                                                  |
| geno_error_coeff    | Slope of genotype error rate. [error] = [offset] + [1-offset]_[coeff]_[1-r2]. Default: 0.0                                                                   |
| r2_info             | INFO field name representing R2 value. Used for representing imputation quality. Default: R2                                                                 |
| min_mac             | Minimum minor allele frequency. Default: 1                                                                                                                   |
| min_callrate        | Minimum call rate. Default: 0.5                                                                                                                              |
| alpha               | Grid of alpha to search for. Default: 0.5                                                                                                                    |
| doublet-prior       | Prior of doublet. Default: 0.5                                                                                                                               |
| demuxlet_out        | Prefix out the demuxlet and dsc-pileup output files. Default: demuxlet_res                                                                                   |

### Genotype-based: Freemuxlet and dsc-pileup

|                            |                                                                                                                                                                |
| -------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| freemuxlet                 | Whether to run Freemuxlet. Default: True                                                                                                                       |
| freemuxlet_preprocess      | Whether to perform pre-processing on the input params.bam for Freemuxlet. True: Perform pre-processing. Otherwise pre-processing is not called. Default: False |
| bam                        | Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed.                                                                                            |
| bai                        | Index of Input SAM/BAM/CRAM file.                                                                                                                              |
| barcodes                   | List of cell barcodes to consider.                                                                                                                             |
| nsample                    | Number of samples multiplexed together                                                                                                                         |
| tag_group                  | Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB Default: CB                  |
| tag_UMI                    | Tag representing UMIs. For 10x genomiucs, use UB. Default: UB                                                                                                  |
| common_variants_freemuxlet | Input VCF/BCF file for dsc-pileup, containing the AC and AN field.                                                                                             |
| sm                         | List of sample IDs to compare to. Default: None (use all)                                                                                                      |
| sm_list                    | File containing the list of sample IDs to compare. Default: None                                                                                               |
| sam_verbose                | Verbose message frequency for SAM/BAM/CRAM. Default: 1000000                                                                                                   |
| vcf_verbose                | Verbose message frequency for VCF/BCF. Default: 10000                                                                                                          |
| skip_umi                   | Do not generate [prefix].umi.gz file, which stores the regions covered by each barcode/UMI pair. Default: False                                                |
| cap_BQ                     | Maximum base quality (higher BQ will be capped). Default: 40                                                                                                   |
| min_BQ                     | Minimum base quality to consider (lower BQ will be skipped). Default: 13                                                                                       |
| min_MQ                     | Minimum mapping quality to consider (lower MQ will be ignored). Default: 20                                                                                    |
| min_TD                     | Minimum distance to the tail (lower will be ignored). Default: 0                                                                                               |
| excl_flag                  | SAM/BAM FLAGs to be excluded. Default: 3844                                                                                                                    |
| min_total                  | Minimum number of total reads for a droplet/cell to be considered. Default: 0                                                                                  |
| min_uniq                   | Minimum number of unique reads (determined by UMI/SNP pair) for a droplet/cell to be considered. Default: 0                                                    |
| min_umi                    | Minimum number of UMIs for a droplet/cell to be considered. Default: 0                                                                                         |
| min_snp                    | Minimum number of SNPs with coverage for a droplet/cell to be considered. Default: 0                                                                           |
| init_cluster               | Input file containing the initial cluster information. Default: None                                                                                           |
| aux_files                  | Turn on writing auxiliary output files. Default: False                                                                                                         |
| verbose                    | Turn on verbose mode with specific verbosity threshold. 0: fully verbose, 100 : no verbose messages. Default: 100                                              |
| doublet_prior              | Prior of doublet. Default: 0.5                                                                                                                                 |
| bf_thres                   | Bayes Factor Threshold used in the initial clustering. Default: 5.41                                                                                           |
| frac_init_clust            | Fraction of droplets to be clustered in the very first round of initial clustering procedure. Default: 0.5                                                     |
| iter_init                  | Iteration for initial cluster assignment (set to zero to skip the iterations). Default: 10                                                                     |
| keep_init_missing          | Keep missing cluster assignment as missing in the initial iteration. Default: False                                                                            |
| freemuxlet_out             | Prefix out the freemuxlet and dsc-pileup output files. Default: freemuxlet_out                                                                                 |

### Genotype-based: Vireo

|                  |                                                                                                                                                                               |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| vireo            | Whether to run Vireo. Default: True                                                                                                                                           |
| vireo_preprocess | Whether to perform pre-processing on the input params.bam for cellSNP-lite. True: Perform pre-processing. Otherwise pre-processing is not called. Default: False              |
| vireo_variant    | Whether to perform cellSNP-lite before running Vireo. True: Run cellSNP-lite. Otherwise cellSNP-lite is not called and params.celldata is used as input. Default: True        |
| celldata         | The cell genotype file in VCF format or cellSNP folder with sparse matrices.                                                                                                  |
| nsample          | Number of donors to demultiplex; can be larger than provided in vcf_donor                                                                                                     |
| vartrixData      | The cell genotype files in vartrix outputs (three/four files, comma separated): alt.mtx,ref.mtx,barcodes.tsv,SNPs.vcf.gz. This will suppress cellData argument. Default: None |
| vcf_donor        | The donor genotype file in VCF format. Default: None                                                                                                                          |
| genoTag          | The tag for donor genotype: GT, GP, PL. Default: GT                                                                                                                           |
| noDoublet        | If use, not checking doublets. Default: False                                                                                                                                 |
| nInit            | Number of random initializations, when GT needs to learn. Default: 50                                                                                                         |
| extraDonor       | Number of extra donor in pre-cluster, when GT needs to learn. Default: 0                                                                                                      |
| extraDonorMode   | Method for searching from extra donors. size: n_cell per donor; distance: GT distance between donor. Default: distance                                                        |
| forceLearnGT     | If use, treat donor GT as prior only. Default: False                                                                                                                          |
| ASEmode          | If use, turn on SNP specific allelic ratio. Default: False                                                                                                                    |
| noPlot           | If use, turn off plotting GT distance. Default: False                                                                                                                         |
| randSeed         | Seed for random initialization. Default: None                                                                                                                                 |
| cellRange        | Range of cells to process, eg. 0-10000. Default: all                                                                                                                          |
| callAmbientRNAs  | If use, detect ambient RNAs in each cell. Default: False                                                                                                                      |
| nproc            | Number of subprocesses for computing, sacrifices memory for speedups. Default: 4                                                                                              |
| vireo_out        | Dirtectory for output files. Default: vireo_out                                                                                                                               |

### Genotype-based: scSplit

|                         |                                                                                                                                                                                                                                                            |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| scSplit                 | Whether to run scSplit. Default: True                                                                                                                                                                                                                      |
| scSplit_preprocess      | Whether to perform pre-processing on the input params.bam for Freebayes and scSplit. True: Perform pre-processing. Otherwise pre-processing is not called. Default: True                                                                                   |
| scSplit_variant         | Whether to perform Freebayes before running scSplit. True: run Freebayes. Otherwise freebayes is not called and params.vcf_mixed is used as input. Default: True                                                                                           |
| vcf_mixed               | VCF from mixed BAM. Default: None                                                                                                                                                                                                                          |
| bam                     | Input Mixed sample BAM.                                                                                                                                                                                                                                    |
| bai                     | Index of mixed sample BAM.                                                                                                                                                                                                                                 |
| barcodes                | Barcodes whitelist.                                                                                                                                                                                                                                        |
| tag_group               | Tag for barcode. Default: CB                                                                                                                                                                                                                               |
| common_variants_scSplit | Common SNVs for scSplit.                                                                                                                                                                                                                                   |
| nsample                 | Expected number of mixed samples.                                                                                                                                                                                                                          |
| refscSplit              | Output Ref count matrix. Default: ref_filtered.csv                                                                                                                                                                                                         |
| altscSplit              | Output Alt count matrix. Default: alt_filtered.csv                                                                                                                                                                                                         |
| subscSplit              | The maximum number of subpopulations in autodetect mode. Default: 10                                                                                                                                                                                       |
| emsscSplit              | Number of EM repeats to avoid local maximum. Default: 30                                                                                                                                                                                                   |
| dblscSplit              | Correction for doublets, Setting to 0 means you would expect no doublets. There will be no refinement on the results if this optional parameter is not specified or specified percentage is less than doublet rates detected during the run. Default: None |
| vcf_donor               | Known individual genotypes to limit distinguishing variants to available variants, so that users do not need to redo genotyping on selected variants.                                                                                                      |
| sample_geno             | Whether to generate sample genotypes based on the split result. Default: True                                                                                                                                                                              |
| scsplit_out             | Dirtectory for scSplit output files. Default: scsplit_out                                                                                                                                                                                                  |

### Genotype-based: Souporcell

|                              |                                                                                                                                                                |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| souporcell                   | Whether to run Souporcell. Default: True                                                                                                                       |
| souporcell_preprocess        | Whether to perform pre-processing on the input params.bam for Souporcell. True: Perform pre-processing. Otherwise pre-processing is not called. Default: False |
| bam                          | Cellranger bam.                                                                                                                                                |
| bai                          | Index of cellranger bam.                                                                                                                                       |
| barcodes                     | Barcodes.tsv from cellranger                                                                                                                                   |
| fasta                        | Reference fasta file.                                                                                                                                          |
| fasta_index                  | Index of reference fasta file.                                                                                                                                 |
| nsample                      | Number of clusters in the BAM file.                                                                                                                            |
| threads                      | Max threads to use. Default: 5                                                                                                                                 |
| ploidy                       | Ploidy, must be 1 or 2. Default: 2                                                                                                                             |
| min_alt                      | Min alt to use locus. Default: 10                                                                                                                              |
| min_ref                      | Min ref to use locus. Default: 10                                                                                                                              |
| max_loci                     | Max loci per cell, affects speed. Default: 2048                                                                                                                |
| restarts                     | Number of restarts in clustering, when there are > 12 clusters we recommend increasing this to avoid local minima. Default: None                               |
| common_variants_souporcell   | Common variant loci or known variant loci vcf, must be vs same reference fasta.                                                                                |
| vcf_donor                    | Known variants per clone in population vcf mode, must be VCF file.                                                                                             |
| known_genotypes_sample_names | Which samples in population vcf from known genotypes option represent the donors in your sample. Default: None                                                 |
| skip_remap                   | Don't remap with minimap2, not recommended unless in conjunction with comman variants. Default: True                                                           |
| ignore                       | Set to True to ignore data error assertions. Default: False                                                                                                    |
| souporcell_out               | Dirtectory for Souporcell output files. Default: souporcell_out                                                                                                |

### Genotype-based: cellSNP-lite

|                         |                                                                                                                                                                                             |
| ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| bam                     | An indexed sam/bam file, comma separated multiple samples.                                                                                                                                  |
| barcodes                | A plain file listing all effective cell barcode.                                                                                                                                            |
| common_variants_cellsnp | A VCF file listing all candidate snps, for fetch each variants.                                                                                                                             |
| targetsVCF              | Similar as regionsVCF, but the next position is accessed by streaming rather than indexing/jumping. Default: None                                                                           |
| sampleList              | A list file containing sample IDs, each per line. Default: None                                                                                                                             |
| sampleIDs               | Comma separated sample ids. Default: None                                                                                                                                                   |
| genotype_cellSNP        | If use, do genotyping in addition to counting. Default: True                                                                                                                                |
| gzip_cellSNP            | If use, the output files will be zipped into BGZF format. Default: True                                                                                                                     |
| printSkipSNPs           | If use, the SNPs skipped when loading VCF will be printed. Default: False                                                                                                                   |
| nproc_cellSNP           | min alt to use locus. Default: 10                                                                                                                                                           |
| refseq_cellSNP          | Faidx indexed reference sequence file. If set, the real (genomic) ref extracted from this file would be used for Mode 2 or for the missing REFs in the input VCF for Mode 1. Default: None. |
| chrom                   | The chromosomes to use, comma separated. Default: None (1-22)                                                                                                                               |
| cellTAG                 | Tag for cell barcodes, turn off with None. Default: CB                                                                                                                                      |
| UMItag                  | Tag for UMI: UB, Auto, None. For Auto mode, use UB if barcodes are inputted, otherwise use None. None mode means no UMI but read counts. Default: Auto                                      |
| minCOUNT                | Minimum aggragated count. Default: 20                                                                                                                                                       |
| minMAF                  | Minimum minor allele frequency. Default: 0.0                                                                                                                                                |
| doubletGL               | If use, keep doublet GT likelihood. Default: False                                                                                                                                          |
| inclFLAG                | Required flags: skip reads with all mask bits unset []. Default: None                                                                                                                       |
| exclFLAG                | Filter flags: skip reads with any mask bits set [UNMAP,SECONDARY,QCFAIL (when use UMI) or UNMAP,SECONDARY,QCFAIL,DUP (otherwise)]. Default: None                                            |
| minLEN                  | Minimum mapped length for read filtering. Default: 30                                                                                                                                       |
| minMAPQ                 | Minimum MAPQ for read filtering. Default: 20                                                                                                                                                |
| maxDEPTH                | Maximum depth for one site of one file (excluding those filtered reads), avoids excessive memory usage; 0 means highest possible value. Default: 0                                          |
| countORPHAN             | If use, do not skip anomalous read pairs. Default: False                                                                                                                                    |
| cellsnp_out             | Dirtectory for cellSNP-lite output files. Default: cellSNP_out                                                                                                                              |

### Genotype-based: Freebayes

|                                 |                                                                                                                                                                                                                                                            |
| ------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| bam                             | Input BAM file to be analyzed.                                                                                                                                                                                                                             |
| bai                             | Index of input BAM file to be analyzed.                                                                                                                                                                                                                    |
| fasta                           | A reference sequence for analysis.                                                                                                                                                                                                                         |
| fasta_index                     | The index of the reference sequence for analysis.                                                                                                                                                                                                          |
| stdin                           | Read BAM input on stdin. Default: False                                                                                                                                                                                                                    |
| targets                         | Limit analysis to targets listed in the BED-format file. Default: None                                                                                                                                                                                     |
| region                          | Limit analysis to the specified chromosome region, 0-base coordinates. If set to None, all chromosomes are considered. Default: None.                                                                                                                      |
| samples                         | Limit analysis to samples listed (one per line) in the file. By default FreeBayes will analyze all samples in its input BAM files. Default: None                                                                                                           |
| populations                     | Each line of FILE should list a sample and a population which it is part of. The population-based bayesian inference model will then be partitioned on the basis of the populations. Default: None                                                         |
| cnv_map                         | Read a copy number map from the BED file. Default: None                                                                                                                                                                                                    |
| vcf_freebayes                   | Name of output VCF file, must be end with .vcf. Default: vcf_freebayes_output.vcf                                                                                                                                                                          |
| gvcf                            | Write gVCF output, which indicates coverage in uncalled regions. Default: False                                                                                                                                                                            |
| gvcf_chunk                      | When writing gVCF output emit a record for every specified number of bases. Default: None                                                                                                                                                                  |
| gvcf_dont_use_chunk             | When writing gVCF output don't emit a record for every specified number of bases. Default: None                                                                                                                                                            |
| variant_input                   | Use variants reported in VCF file as input to the algorithm. Variants in this file will included in the output even if there is not enough support in the data to pass input filters. Default: None                                                        |
| only_use_input_alleles          | Only provide variant calls and genotype likelihoods for sites and alleles which are provided in the VCF input, and provide output in the VCF for all input alleles, not just those which have support in the data. Default: False                          |
| haplotype_basis_alleles         | When specified, only variant alleles provided in this input VCF will be used for the construction of complex or haplotype alleles. Default: None                                                                                                           |
| report_all_haplotype_alleles    | At sites where genotypes are made over haplotype alleles, provide information about all alleles in output, not only those which are called. Default: False                                                                                                 |
| report_monomorphic              | Report even loci which appear to be monomorphic, and report all considered alleles, even those which are not in called genotypes. Default: False                                                                                                           |
| pvar                            | Report sites if the probability that there is a polymorphism at the site is greater than N. Default: 0.0                                                                                                                                                   |
| strict_vcf                      | Generate strict VCF format (FORMAT/GQ will be an int). Default: False                                                                                                                                                                                      |
| theta                           | The expected mutation rate or pairwise nucleotide diversity among the population under analysis. This serves as the single parameter to the Ewens Sampling Formula prior model. Default: 0.001                                                             |
| ploidy                          | Sets the default ploidy for the analysis. Default: 2                                                                                                                                                                                                       |
| pooled_discrete                 | Assume that samples result from pooled sequencing. Model pooled samples using discrete genotypes across pools. When using this flag, set --ploidy to the number of alleles in each sample or use the --cnv-map to define per-sample ploidy. Default: False |
| pooled_continuous               | Output all alleles which pass input filters, regardless of genotyping outcome or model. Default: False                                                                                                                                                     |
| use_reference_allele            | This flag includes the reference allele in the analysis as if it is another sample from the same population. Default: False                                                                                                                                |
| reference_quality               | Assign mapping quality to the reference allele at each site and base quality. Default: 100,60                                                                                                                                                              |
| no_snps                         | Ignore SNP alleles. Default: False                                                                                                                                                                                                                         |
| no_indels                       | Ignore insertion and deletion alleles. Default: True                                                                                                                                                                                                       |
| no_mnps                         | Ignore multi-nuceotide polymorphisms, MNPs. Default: True                                                                                                                                                                                                  |
| no_complex                      | Ignore complex events (composites of other classes). Default: True                                                                                                                                                                                         |
| use_best_n_alleles              | Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores. Set to 0 to use all. Default: 0                                                                                                                                          |
| haplotype_length                | Allow haplotype calls with contiguous embedded matches of up to this length. Set N=-1 to disable clumping. Default: 3                                                                                                                                      |
| min_repeat_size                 | When assembling observations across repeats, require the total repeat length at least this many bp. Default: 5                                                                                                                                             |
| min_repeat_entropy              | To detect interrupted repeats, build across sequence until it has entropy > N bits per bp. Set to 0 to turn off. Default: 1                                                                                                                                |
| no_partial_observations         | Exclude observations which do not fully span the dynamically-determined detection window. Default: None, to use all observations, dividing partial support across matching haplotypes when generating haplotypes.                                          |
| dont_left_align_indels          | Turn off left-alignment of indels, which is enabled by default. Default: False                                                                                                                                                                             |
| use_duplicate_reads             | Include duplicate-marked alignments in the analysis. Default: False, to exclude duplicates marked as such in alignments                                                                                                                                    |
| min_mapping_quality             | Exclude alignments from analysis if they have a mapping quality less than Q. Default: 1                                                                                                                                                                    |
| min_base_quality                | Exclude alleles from analysis if their supporting base quality is less than Q. Default: 1                                                                                                                                                                  |
| min_supporting_allele_qsum      | Consider any allele in which the sum of qualities of supporting observations is at least Q. Default: 0                                                                                                                                                     |
| min_supporting_mapping_qsum     | Consider any allele in which and the sum of mapping qualities of supporting reads is at least Q. Default: 0                                                                                                                                                |
| mismatch_base_quality_threshold | Count mismatches toward --read-mismatch-limit if the base quality of the mismatch is >= Q. Default: 10                                                                                                                                                     |
| read_mismatch_limit             | Exclude reads with more than N mismatches where each mismatch has base quality >= mismatch-base-quality-threshold. Default: None, ~unbounded                                                                                                               |
| read_max_mismatch_fraction      | Exclude reads with more than N [0,1] fraction of mismatches where each mismatch has base quality >= mismatch-base-quality-threshold. Default: 1.0                                                                                                          |
| read_snp_limit                  | Exclude reads with more than N base mismatches, ignoring gaps with quality >= mismatch-base-quality-threshold. Default: None, ~unbounded                                                                                                                   |
| read_indel_limit                | Exclude reads with more than N separate gaps. Default: None, ~unbounded                                                                                                                                                                                    |
| standard_filters                | Use stringent input base and mapping quality filters equivalent to -m 30 -q 20 -R 0 -S 0. Default: False                                                                                                                                                   |
| min_alternate_fraction          | Require at least this fraction of observations supporting an alternate allele within a single individual in in order to evaluate the position. Default: 0.05                                                                                               |
| min_alternate_count             | Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position. Default: 2                                                                                                        |
| min_alternate_qsum              | Require at least this sum of quality of observations supporting an alternate allele within a single individual in order to evaluate the position. Default: 0                                                                                               |
| min_alternate_total             | Require at least this count of observations supporting an alternate allele within the total population in order to use the allele in analysis. Default: 1                                                                                                  |
| min_coverage                    | Require at least this coverage to process a site. Default: 0                                                                                                                                                                                               |
| max_coverage                    | Do not process sites with greater than this coverage. Default: None, no limit                                                                                                                                                                              |
| no_population_priors            | Equivalent to --pooled-discrete --hwe-priors-off and removal of Ewens Sampling Formula component of priors. Default: False                                                                                                                                 |
| hwe_priors_off                  | Disable estimation of the probability of the combination arising under HWE given the allele frequency as estimated by observation frequency. Default: False                                                                                                |
| binomial_obs_priors_off         | Disable incorporation of prior expectations about observations. Uses read placement probability, strand balance probability, and read position (5'-3') probability. Default: False                                                                         |
| allele_balance_priors_off       | Disable use of aggregate probability of observation balance between alleles as a component of the priors. Default: False                                                                                                                                   |
| observation_bias                | Read length-dependent allele observation biases from the file. Default: None                                                                                                                                                                               |
| base_quality_cap                | Limit estimated observation quality by capping base quality at Q. Default: None                                                                                                                                                                            |
| prob_contamination              | An estimate of contamination to use for all samples. Default: 10e-9                                                                                                                                                                                        |
| legacy_gls                      | Use legacy (polybayes equivalent) genotype likelihood calculations. Default: False                                                                                                                                                                         |
| contamination_estimates         | A file containing per-sample estimates of contamination, such as those generated by VerifyBamID. Default: None                                                                                                                                             |
| report_genotype_likelihood_max  | Report genotypes using the maximum-likelihood estimate provided from genotype likelihoods. Default: False                                                                                                                                                  |
| genotyping_max_iterations       | Iterate no more than N times during genotyping step. Default: 1000                                                                                                                                                                                         |
| genotyping_max_banddepth        | Integrate no deeper than the Nth best genotype by likelihood when genotyping. Default: 6                                                                                                                                                                   |
| posterior_integration_limits    | Integrate all genotype combinations in our posterior space which include no more than N samples with their Mth best data likelihood. Default: 1,3                                                                                                          |
| exclude_unobserved_genotypes    | Skip sample genotypings for which the sample has no supporting reads. Default: False                                                                                                                                                                       |
| genotype_variant_threshold      | Limit posterior integration to samples where the second-best genotype likelihood is no more than log(N) from the highest genotype likelihood for the sample. Default: None, ~unbounded                                                                     |
| use_mapping_quality             | Use mapping quality of alleles when calculating data likelihoods. Default: False                                                                                                                                                                           |
| harmonic_indel_quality          | Use a weighted sum of base qualities around an indel, scaled by the distance from the indel. Default: False, use a minimum BQ in flanking sequence.                                                                                                        |
| read_dependence_factor          | Incorporate non-independence of reads by scaling successive observations by this factor during data likelihood calculations. Default: 0.9                                                                                                                  |
| genotype_qualities              | Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output Default: False                                                                                                                                     |
| debug                           | Print debugging output. Default: False                                                                                                                                                                                                                     |
| dd                              | Print more verbose debugging output (requires "make DEBUG"). Default: False                                                                                                                                                                                |

### Donor matching

|                       |                                                                                                                                                                                                                                                                         |
| --------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| match_donor           | Whether to match donors. Default: True                                                                                                                                                                                                                                  |
| demultiplexing_result | A CSV file with demultiplexing assignment when running in donor_match mode. In other modes, the input is passed by the pipeline automatically. Default: None                                                                                                            |
| match_donor_method1   | The method name to match donors. If None, all genotype-based methods are compared. Default: None                                                                                                                                                                        |
| match_donor_method2   | The method name to match donors. If None, all hashing-based methods are compared. Default: None                                                                                                                                                                         |
| findVariants          | Whether to extract a subset of informative variants when best genotype-based method for donor matching is vireo. `default`: subset as described in paper; `vireo`: subset by Vireo; `True`: subset using both methods; `False`: not extracting variants. Default: False |
| variant_count         | The threshold for the minimal read depth of a variant in the cell group when subseting the informative variants by default. Default: 10                                                                                                                                 |
| variant_pct           | The threshold for the minimal frequency of the alternative or reference allele to determine the dominant allele of a variant in the cell group when subseting the informative variants by default. Default: 0.9                                                         |
| vireo_parent_dir      | A parent folder which contains the output folder of vireo in the format of `vireo_[taskID/sampleId]` generated by hadge pipeline when running in donor_match mode. In other modes, the input is passed by the pipeline automatically. Default: None                     |
