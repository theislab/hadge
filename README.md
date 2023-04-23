# hagen

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

![Caption](pipeline_v2.png)
## **Introduction**
Hagen is a one-stop pipeline for demultiplexing single cell mixtures. It consists of 14 methods across two workflows: hashing and genetic demultiplexing, which can be run in 3 modes. 

The genetic demultiplexing workflow includes 5 methods: 
* Freemuxlet 
* Demuxlet
* Vireo
* Souporcell
* scSplit

The hashing workflow includes 9 methods: 
* hashedDrops
* Multiseq
* HTODemux
* Demuxem
* Solo
* HashSolo
* TODO

## **Installation**

The Hagen pipeline is implemented in Nextflow. To get started, you need to install Nextflow. Please refer to [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) for more details. Alternatively, you can also install Nextflow via [conda](https://anaconda.org/bioconda/nextflow).

As next, please clone the repository
```
git clone https://github.com/theislab/hagen.git
cd hagen
```

If you want to use Souporcell later, you should also download the singularity image in the project directory with the command `singularity pull shub://wheaton5/souporcell`.


To run the pipeline, you need to prepare the input data and set the configuration file correctly.

## **Data preparation**
The input data depends heavily on the demultiplexing tools. In the following table, you will find the minimal input data required by different tools.

### Genetic demultiplexing workflow

| Demultiplexing method 	| Input data                                                                           	|
|-----------------------	|--------------------------------------------------------------------------------------	|
| Demuxlet              	| - Alignment (BAM)<br>- Barcode (TSV)<br>- Genotype reference (VCF)                   	|
| Freemuxlet            	| - Alignment (BAM)<br>- Barcode (TSV)<br>- Genotypes from referenced population (VCF) 	|
| Vireo                 	| - Cell genotype (VCF or cellSNP folder)                                              	|
| Souporcell            	| - Alignment (BAM)<br>- Barcode (TSV)<br>- Reference genome (FASTA)                   	|
| scSplit               	| - Alignment (BAM)<br>- Barcode (TSV)<br>- Genotypes from referenced population (VCF) 	|
|                       	|                                                                                      	|

You may see that some tools share some input data in common, so we have only one parameter for the same input.

| Input data                                 	| Parameter                              	|
|--------------------------------------------	|----------------------------------------	|
| Alignment (BAM)                            	| `params.bam`<br>`params.bai`           	|
| Barcode (TSV)                              	| `params.barcodes`                      	|
| Genotype reference (VCF)                   	| `params.vcf_donor`                     	|
| Genotypes from referenced population (VCF) 	| `params.vcf_mixed`                     	|
| Reference genome (FASTA)                   	| `params.fasta`<br>`params.fasta_index` 	|
| Cell genotype (VCF or cellSNP folder)      	| `params.celldata`                      	|

In case you want to perform genetic demultiplexing on pre-processed data, we provide a process for this purpose, which is in concordance with [the instruction of scSplit](https://github.com/jon-xu/scSplit). It only requires the Alignment (BAM) file as input. To activate the pre-processing step, set `[method]_preprocess` as `True`.

In case you don't have any cell genotypes or variants called from mixed samples yet, we provide two processes for variant calling. 
| Variant calling methods 	| Input data                                                  	| Parameter                                                                	| Output                      	|
|-------------------------	|-------------------------------------------------------------	|--------------------------------------------------------------------------	|-----------------------------	|
| freebayes               	| - Alignment (BAM)<br>- Reference genome (FASTA)             	| `params.bam`<br>`params.bai`<br>`params.fasta`<br>`params.fasta_index`   	| Variants from mixed samples 	|
| cellsnp-lite            	| - Alignment (BAM)<br>- Barcode (TSV)<br>- Common SNPs (VCF) 	| `params.bam`<br>`params.bai`<br>`params.barcodes`<br>`params.regionsVCF` 	| Cell genotypes              	|

For benchmarking, you can have following options for `scsplit_variant`:
* `False`: inactivate variant calling, get the input data from `params.vcf_mixed`
* `freebayes`: activate freebayes
* Otherwise: activate freebayes, take the output of freebayes and use `params.vcf_mixed` as well

For benchmarking, you can have following options for `scsplit_variant`:
* `False`: inactivate variant calling, get the input data from `params.celldata`
* `cellsnp`: activate cellsnp
* Otherwise: activate cellsnp, take the output of freebayes and use `params.celldata` as well

### Hashing workflow

| Demultiplexing method 	| Input data                                                                                                         	| Parameter                                                  	|
|-----------------------	|--------------------------------------------------------------------------------------------------------------------	|------------------------------------------------------------	|
| HTODemux              	| - Seurat object with both UMI and hashing count matrix (RDS)                                                       	| `params.rdsObj_htodemux`                                   	|
| Multiseq              	| - Seurat object with both UMI and hashing count matrix (RDS)                                                       	| `params.rdsObj_htodemux`                                   	|
| Solo                  	| - 10x mtx directory with UMI count matrix (Directory)                                                              	| `params.rna_matrix_solo`                                   	|
| HashSolo              	| - HDF5 file with hashing count matrix (H5)                                                                         	| `params.hto_h5_hashsolo`                                   	|
| HashedDrops           	| - 10x mtx directory with hashing count matrix (Directory)                                                          	| `params.hto_matrix_hashedDrops`                            	|
| Demuxem               	| - 10x mtx directory with UMI count matrix (Directory)<br>- 10x mtx directory with hashing count matrix (Directory) 	| `params.hto_matrix_demuxem`<br>`params.rna_matrix_demuxem` 	|

Similar as in the genetic demultiplexing workflow, we provide a pre-processing step specifically for HTODemux and Multiseq. The input can be either an RDS object via `params.rdsObject` or the UMI and hashing count matrix `params.umi_matrix_preprocess` and `params.hto_matrix_preprocess`.

For benchmarking, you can have following options for `[htodemux/multiseq]_preprocess`:
* `True`: activate pre-proecessing
* `False`: inactivate pre-proecessing, get the input data from `params.rdsObj_[method]`
* Otherwise: activate pre-proecessing, take the output and use `params.rdsObj_[method]` as well

  
## **Pipeline configuration**
In the `nextflow.config` file, 
* Set the `mode`: genetic, hashing or rescue.
* Define the folder name `outdir` to save the output files. This will create a folder automatically in the project directory. The result of each process will be saved in the folder `$projectDir/$params.outdir/$params.mode`

* We strongly recommend to run each process in a separate container or in a conda environment. For example, you have three options to set a conda envrionment for each progress:
    ```
    // dont forget to enable conda
    conda.enable = true
    process {
        // Use Conda environment files
        withName:scSplit {
            conda = './conda/scsplit.yml' 
        }
        // Use Conda package names
        withName:cellSNP {
            conda = 'bioconda::cellsnp-lite'
        }
        // Use existing Conda environments
        withName:summary {
            conda = '/path/to/an/existing/env/directory'
        }
    }

    ```
* The pipeline can be run either locally or on a HPC. You can set the executor by running the pipeline with `-profile standard` or `-profile cluster`. Of course, you can add other profiles if you want.
    ```
    profiles{
        standard {
            process.executor = 'local'
        }

        cluster {
            process.executor = 'slurm'
        }
    }
    ```
* Feel free to add other configurations, e.g. the number of CPUS, the memory allocation, etc. If you are new to Nextflow framework, please visit the [Nextlfow page](https://www.nextflow.io/docs/latest/config.html#).

Finally, you can run the pipeline with: 

    nextflow run main.nf -profile standard

## Pipeline output

### Intermediate Output
Each demultiplexing process will generate some intermediate files in the folder in the format `[method]/[method]_[task_ID]`, e.g. `htodemux/htodemux_1`. In this folder, you can find following files:
* `params.csv`: specified parameters in the task
* Output of the task (depends on the tool)

### Final output
After each demultiplexing workflow, the pipeline will generate some TSV files to summarize the results in the folder `[workflow]_summary`.
*  `[method]_classification.csv`: classification of all trials for a given method
    |  Barcode  	| multiseq_1 	| multiseq_2 	| ... 	|
    |:---------:	|:----------:	|:----------:	|:---:	|
    | barcode-1 	|   singlet  	|   singlet  	| ... 	|
    | barcode-2 	|   doublet  	|  negative  	| ... 	|
    |    ...    	|     ...    	|     ...    	| ... 	|
*  `[method]_assignment.csv`: assignment of all trials for a given method
    |  Barcode  	| multiseq_1 	| multiseq_2 	| ... 	|
    |:---------:	|:----------:	|:----------:	|:---:	|
    | barcode-1 	|   donor-1  	|   donor-2  	| ... 	|
    | barcode-2 	|   doublet  	|  negative  	| ... 	|
    |    ...    	|     ...    	|     ...    	| ... 	|
* `[method]_params.csv`: specified paramters of all trials for a given method
    |    Argument   	|     Value     |
    |    :---------:	|  :----------:	|
    |  seuratObejctPath |      Path 	|
    |     quantile   	|      0.7  	|
    |       ...     	|      ...    	|
* `[workflow]_classification_all.csv`: classification of all trials across different methods
    |  Barcode  	| multiseq_1 	| htodemux_1 	| ... 	|
    |:---------:	|:----------:	|:----------:	|:---:	|
    |    ...    	|     ...    	|     ...    	| ... 	|
*  `[workflow]_assignment_all.csv`: save the assignment of all trials across different methods
    |  Barcode  	| multiseq_1 	| htodemux_1 	| ... 	|
    |:---------:	|:----------:	|:----------:	|:---:	|
    |    ...    	|     ...    	|     ...    	| ... 	|
In the `rescue` mode: 
* The pipeline merges the results of hashing and genetic demultiplexing tools into `classification_all_genetic_and_hash.csv` and `assignment_all_genetic_and_hash.csv` in the `summary` folder.
* Donor matching: in the folder `donor_match/donor_match_[method1]_[method2]` you will find:
    * folder`[method1]_[trial_ID]_vs_[method2]_[trial_ID]` with:
        * `correlation_res.csv`: correlation scores of donor matching
        `concordance_heatmap.png`: a heatmap visualising the the correlation scores
        *  `donor_match.csv`: A map between hashtag and donor identity.
    * Once the optimal match among all trials between `method1` and `method2` is found, the pipeline maps the result of `method1` to `method2` and provides:
        * `all_assignment_after_match.csv`: Assignment of all cell barcodes
        * `intersect_assignment_after_match.csv`: Assignment of joint cell barcodes
    * If `method1` of the optimal match is `vireo` and identification of donor-specific variants is enabled:
        *  `donor_specific_variants.csv`: a list of donor-specific variants
        * `donor_specific_variants_upset.png`: An upset plot showing the number of donor-specific variants
        * `donor_genotype_subset_by_default_matched.vcf`: Donor genotypes of donor-specific variants



## Parameters

## Credits

## Citations