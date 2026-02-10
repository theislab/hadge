# nf-core/hadge: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/hadge/usage](https://nf-co.re/hadge/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

### The rescue mode

The joint call of hashing and genetic deconvolution methods has been shown to be beneficial for cell recovery rate and calling accuracy.
hadge provides a rescue mode to run both genotype- and hashing-based approaches jointly to rescue problematic hashing experiments in cases where donors are genetically distinct.
In this scenario, samples of both hashing and genetic multiplexing experiments are deconvoluted simultaneously.
Furthermore, hadge allows for the automatic determination of the best combination of hashing and SNP-based donor deconvolution tools.

```csv title="samplesheet.csv"
sample,rna_matrix,hto_matrix,bam,vcf,n_samples,barcodes
id1,rna.tar.gz,hto.tar.gz,chr21.bam,donor_genotype_chr21.vcf,2,barcodes.tsv
id2,rna.tar.gz,hto.tar.gz,chr21.bam,donor_genotype_chr21.vcf,2,barcodes.tsv
id3,rna.tar.gz,hto.tar.gz,chr21.bam,donor_genotype_chr21.vcf,2,barcodes.tsv
```

Now, you can run the pipeline using:

```bash
nextflow run nf-core/hadge \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --mode rescue \
   --hash_tools htodemux,hasheddrops,multiseq,gmm-demux,bff,hashsolo \
   --genetic_tools demuxlet,freemuxlet,vireo,souporcell \
   --fasta <FASTADIR>
```

### The genetic mode

Genotyped-based deconvolution assigns cells to donors using genetic variation.
This can be performed with donor genotypes or, if these are unavailable, using reference panels in genotype-free mode (e.g., 1000 Genomes).
Finally, it assigns SNPs to cells to determine donor identity but requires additional genotyping.

```csv title="samplesheet.csv"
sample,bam,vcf,n_samples,barcodes
id1,donor_genotype_chr21.vcf,2,barcodes.tsv
id2,donor_genotype_chr21.vcf,2,barcodes.tsv
id3,chr21.bam,donor_genotype_chr21.vcf,2,barcodes.tsv
```

Now, you can run the pipeline using:

```bash
nextflow run nf-core/hadge \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --mode genetic \
   --genetic_tools demuxlet,freemuxlet,vireo,souporcell \
   --fasta <FASTADIR>
```

:::info
A FASTA file is only required if `--genetic_tools` includes souporcell.
If a FASTA file is unavailable, you can specify the organism using `--genome`, and the pipeline will download the full reference genome automatically.
However, to avoid long download times and high bandwidth usage, we recommend providing your own local reference genome with `--fasta`.
:::

### The hashing mode

Cell hashing tags cells with unique oligo barcodes so samples can be pooled.
Separate scRNA and HTO libraries are sequenced, producing count matrices used to determine each cell’s sample of origin.

```csv title="samplesheet.csv"
sample,rna_matrix,hto_matrix
id1,rna.tar.gz,hto.tar.gz
id2,rna.tar.gz,hto.tar.gz
id3,rna.tar.gz,hto.tar.gz
```

Now, you can run the pipeline using:

```bash
nextflow run nf-core/hadge \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --mode hashing
   --hash_tools htodemux,hasheddrops,multiseq,gmm-demux,bff,hashsolo \
```

### The donor match mode

This mode utilizes the donor matching component from the rescue mode, but requires manual input for several stages.
To run all steps of donor matching, you must provide the demultiplexing results, filtered variants, and both cell and donor genotypes.
For detailed specifications on these input parameters, refer to the [parameter documentation](https://nf-co.re/hadge/parameters).

```csv title="samplesheet.csv"
sample,n_samples
id1,2
```

Now, you can run the pipeline using:

```bash
nextflow run nf-core/hadge \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --mode donor_match
   --demultiplexing_result <DIR> \
   --vireo_filtered_variants <DIR> \
   --cell_genotype <DIR> \
   --gt_donors <DIR> \
```

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline.
Use this parameter to specify its location.
It has to be a comma-separated file with a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

Each row in the sample sheet represents a distinct single-cell multiplexing experiment.
The `sample` column must contain a unique identifier for each experiment.
This format allows you to process multiple deconvolutions in a single run.
While a full example is provided below, some columns may be optional depending on the mode you select.

```csv title="samplesheet.csv"
sample,rna_matrix,hto_matrix,bam,vcf,n_samples,barcodes
id1,rna.tar.gz,hto.tar.gz,chr21.bam,donor_genotype_chr21.vcf,2,barcodes.tsv
id2,rna.tar.gz,hto.tar.gz,chr21.bam,donor_genotype_chr21.vcf,2,barcodes.tsv
id3,rna.tar.gz,hto.tar.gz,chr21.bam,donor_genotype_chr21.vcf,2,barcodes.tsv
```

| Column       | Description                                                                                                                                                                            |
| ------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`     | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `rna_matrix` | Full path to the RNA-Seq count matrices provided in a 10x Genomics format and compressed as `.tar.gz`.                                                                                 |
| `hto_matrix` | Full path to the hashing count matrices provided in a 10x Genomics format and compressed as `.tar.gz`.                                                                                 |
| `bam`        | Full path to the alignment file (`.bam`).                                                                                                                                              |
| `vcf`        | Full path to the list of common SNPs (`.vcf`).                                                                                                                                         |
| `n_samples`  | The number of multiplexed donors.                                                                                                                                                      |
| `barcodes`   | Full path to the list of cell barcodes (e.g., `barcodes.tsv` from Cell Ranger)                                                                                                         |

:::tip{collapse title="Samplesheet Input Requirements by Module"}

| Mode        | sample | rna_matrix | hto_matrix | bam | barcodes | n_samples | vcf |
| ----------- | :----: | :--------: | :--------: | :-: | :------: | :-------: | :-: |
| rescue      |   ✅   |     ✅     |     ✅     | ✅  |    ✅    |    ✅     | ✅  |
| genetic     |   ✅   |     ❌     |     ❌     | ✅  |    ✅    |    ✅     | ✅  |
| hashing     |   ✅   |     ✅     |     ✅     | ❌  |    ❌    |    ❌     | ❌  |
| donor_match |   ✅   |     ❌     |     ❌     | ❌  |    ❌    |    ✅     | ❌  |

| Module      | sample | rna_matrix | hto_matrix | bam | barcodes | n_samples | vcf |
| ----------- | :----: | :--------: | :--------: | :-: | :------: | :-------: | :-: |
| htodemux    |   ✅   |     ✅     |     ✅     | ❌  |    ❌    |    ❌     | ❌  |
| multiseq    |   ✅   |     ✅     |     ✅     | ❌  |    ❌    |    ❌     | ❌  |
| bff         |   ✅   |     ❌     |     ✅     | ❌  |    ❌    |    ❌     | ❌  |
| demuxem     |   ✅   |     ✅     |     ✅     | ❌  |    ❌    |    ❌     | ❌  |
| gmm-demux   |   ✅   |     ❌     |     ✅     | ❌  |    ❌    |    ❌     | ❌  |
| hasheddrops |   ✅   |    ✅\*    |     ✅     | ❌  |    ❌    |    ❌     | ❌  |
| hashsolo    |   ✅   |     ❌     |     ✅     | ❌  |    ❌    |    ❌     | ❌  |
| vireo       |   ✅   |     ❌     |     ❌     | ✅  |    ✅    |    ✅     | ✅  |
| demuxlet    |   ✅   |     ❌     |     ❌     | ✅  |    ❌    |    ❌     | ✅  |
| freemuxlet  |   ✅   |     ❌     |     ❌     | ✅  |    ❌    |    ✅     | ✅  |
| souporcell  |   ✅   |     ❌     |     ❌     | ✅  |    ✅    |    ✅     | ❌  |

\* if `params.hasheddrops_runEmptyDrops` is true
:::

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/hadge --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/hadge -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/hadge
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/hadge releases page](https://github.com/nf-core/hadge/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
