# Usage
## **Input data preparation**
The input data depends heavily on the deconvolution tools. In the following table, you will find the minimal input data required by different tools.

### Genetics-based deconvolution methods:

| Deconvolution methods 	| Input data                                                                           	|
|-----------------------	|--------------------------------------------------------------------------------------	|
| Demuxlet              	| - Alignment (BAM)<br>- Barcode (TSV)<br>- Genotype reference (VCF)                   	|
| Freemuxlet            	| - Alignment (BAM)<br>- Barcode (TSV)<br>- Genotypes from referenced population (VCF) 	|
| Vireo                 	| - Cell genotype (VCF or cellSNP folder)                                              	|
| Souporcell            	| - Alignment (BAM)<br>- Barcode (TSV)<br>- Reference genome (FASTA)                   	|
| scSplit               	| - Alignment (BAM)<br>- Barcode (TSV)<br>- Genotypes from referenced population (VCF) 	|
|                       	|                                                                                      	|

You may see that some tools share some input data in common, so we set only one parameter for the same input for benchmarking.

| Input data                                 	| Parameter                              	|
|--------------------------------------------	|----------------------------------------	|
| Alignment (BAM)                            	| `params.bam`<br>`params.bai`           	|
| Barcode (TSV)                              	| `params.barcodes`                      	|
| Genotype reference (VCF)                   	| `params.vcf_donor`                     	|
| Genotypes from referenced population (VCF) 	| `params.vcf_mixed`                     	|
| Reference genome (FASTA)                   	| `params.fasta`<br>`params.fasta_index` 	|
| Cell genotype (VCF or cellSNP folder)      	| `params.celldata`                      	|


#### <strong> Pre-processing </strong>
In case you want to perform genetics-based deconvolution on pre-processed data, we provide a process in concordance with [the instruction of scSplit](https://github.com/jon-xu/scSplit). It only requires the Alignment (BAM) file as input. To specify which method is performed on the pre-processed data : set `[method]_preprocess = True`.


#### <strong> Variant calling </strong> 

In case you don't have any cell genotypes or variants called from mixed samples yet, we provide two processes for variant calling. 
| Variant calling methods 	| Input data                                                  	| Parameter                                                                	| Output                      	|
|-------------------------	|-------------------------------------------------------------	|--------------------------------------------------------------------------	|-----------------------------	|
| freebayes               	| - Alignment (BAM)<br>- Reference genome (FASTA)             	| `params.bam`<br>`params.bai`<br>`params.fasta`<br>`params.fasta_index`   	| Variants from mixed samples 	|
| cellsnp-lite            	| - Alignment (BAM)<br>- Barcode (TSV)<br>- Common SNPs (VCF) 	| `params.bam`<br>`params.bai`<br>`params.barcodes`<br>`params.regionsVCF` 	| Cell genotypes              	|

You can have following options for `scsplit_variant`:
* `False`: inactivate variant calling, get the input data from `params.vcf_mixed`
* `freebayes`: activate freebayes
* Otherwise: activate freebayes, take the output of freebayes and use `params.vcf_mixed` as well

You can have following options for `scsplit_variant`:
* `False`: inactivate variant calling, get the input data from `params.celldata`
* `cellsnp`: activate cellsnp
* Otherwise: activate cellsnp, take the output of freebayes and use `params.celldata` as well


### Hashing-based deconvolution workflow

| Deconvolution method 	| Input data                                                                                                         	| Parameter                                                  	|
|-----------------------	|--------------------------------------------------------------------------------------------------------------------	|------------------------------------------------------------	|
| HTODemux              	| - Seurat object with both UMI and hashing count matrix (RDS)                                                       	| `params.rdsObj_htodemux`                                   	|
| Multiseq              	| - Seurat object with both UMI and hashing count matrix (RDS)                                                       	| `params.rdsObj_htodemux`                                   	|
| Solo                  	| - 10x mtx directory with UMI count matrix (Directory)                                                              	| `params.rna_matrix_solo`                                   	|
| HashSolo              	| - HDF5 file with hashing count matrix (H5)                                                                         	| `params.hto_h5_hashsolo`                                   	|
| HashedDrops           	| - 10x mtx directory with hashing count matrix (Directory)                                                          	| `params.hto_matrix_hashedDrops`                            	|
| Demuxem               	| - 10x mtx directory with UMI count matrix (Directory)<br>- 10x mtx directory with hashing count matrix (Directory) 	| `params.hto_matrix_demuxem`<br>`params.rna_matrix_demuxem` 	|


#### <strong> Pre-processing </strong>

Similar as in the genetic demultiplexing workflow, we provide a pre-processing step specifically for HTODemux and Multiseq. The input can be either an RDS object via `params.rdsObject` or the UMI and hashing count matrix `params.umi_matrix_preprocess` and `params.hto_matrix_preprocess`.

For benchmarking, you can have following options for `[htodemux/multiseq]_preprocess`:
* `True`: activate pre-proecessing
* `False`: inactivate pre-proecessing, get the input data from `params.rdsObj_[method]`
* Otherwise: activate pre-proecessing, take the output and use `params.rdsObj_[method]` as well


## **Pipeline configuration**
### **Conda environments:**
We provide a  `environment.yml` file for each process. But you can also use local conda environments to run a process:

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

### **Containers:**
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
// dont forget to enable docker

docker {
    enabled = true
}

```


### **Executor and resource specifications:**
* The pipeline can be run either locally or on a HPC. You can set the executor by running the pipeline with `-profile standard` or `-profile cluster`. Of course, you can add other profiles if you want. 
* Feel free to add other configurations, e.g. the number of CPUS, the memory allocation, etc. If you are new to Nextflow framework, please visit the [Nextlfow page](https://www.nextflow.io/docs/latest/config.html#).
*  As default, the pipeline is run locally with the standard profile, where all processes annotated with the big_mem label are assigned 4 cpus and 16 Gb of memory.

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
            executor = 'local'
            withLabel: big_mem {
                cpus = 32
                memory = 64.GB
            }
            withLabel: small_mem {
                cpus = 16
                memory = 32.GB
                // queue = ...
            }
        }
    }
}

```

## **Parameters**