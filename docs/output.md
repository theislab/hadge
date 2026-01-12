# nf-core/hadge: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

### Genetic-based deconvolution

This subworkflow runs when the `--mode` is set to `genetic` or `rescue`. It saves the results of each genetic-based deconvolution specified with `--genetic_tools` in a folder named after the corresponding tool. Additionally, this step summarizes all deconvolution results in a single summary table. The assignments CSV provides, for each cell barcode, one of the following categories:

| Category      | Description                                    |
| ------------- | ---------------------------------------------- |
| Cluster label | A cluster label (e.g., `0`, `1`, …)            |
| doublet       | More than one cluster is assigned to this cell |
| negative      | All other cases (e.g., undetected cells)       |

Classification works similarly, except that cluster labels are replaced by the label `singlet`.

<details markdown="1">
<summary>Output files</summary>

- `genetic/`
  - TODO add every tool ...
  - `*.fastp.html`: Trimming report in html format.
  - `*.fastp.json`: Trimming report in json format.
- `genetic/summary/`
  - `/*/*_genetic_summary_(assignment|classification).csv`: Summary of all assigned/classified cells from each genetic-based deconvolution tool, merged into a single table.
  - `/*/*_genetic_overview_(assignment|classification).csv`: This table summarizes each genetic-based deconvolution tool (before merging) by reporting its total barcode count, the number of barcodes it shares with every other method, and the counts of each donor label or classification category (e.g., `0`, `1`, `singlet`, `doublet`, `negative`).

### Hashing-based deconvolution

This subworkflow runs when the `--mode` is set to `hashing` or `rescue`. It saves the results of each hashing-based deconvolution specified with `--hash_tools` in a folder named after the corresponding tool. Additionally, this step summarizes all deconvolution results in a single summary table. The assignments CSV provides, for each cell barcode, one of the following categories:

| Category    | Description                                                       |
| ----------- | ----------------------------------------------------------------- |
| Donor label | A HTO label that identifies the donor (e.g., `HTO-1`, `HTO-2`, …) |
| doublet     | More than one cluster is assigned to this cell                    |
| negative    | All other cases (e.g., undetected cells)                          |

Classification works similarly, except that cluster labels are replaced by the label `singlet`.

<details markdown="1">
<summary>Output files</summary>

- `hashing/`
  - TODO add every tool ...
  - `*.fastp.html`: Trimming report in html format.
  - `*.fastp.json`: Trimming report in json format.
- `hashing/summary/`
  - `/*/*_hashing_summary_(assignment|classification).csv`: Summary of all assigned/classified cells from each hashing-based deconvolution tool, merged into a single table.
  - `/*/*_hashing_overview_(assignment|classification).csv`: This table summarizes each hashing-based deconvolution tool (before merging) by reporting its total barcode count, the number of barcodes it shares with every other method, and the counts of each donor label or classification category (e.g., `HTO-1`, `HTO-2`, `singlet`, `doublet`, `negative`).

</details>

### Donor matching

<!-- TODO add output files -->

### Summary

When running the pipeline in `rescue` mode, the combined assignment/classification table from all used genetic- and hashing-based deconvolution tools will be generated here. Additionally, an AnnData/MuData object is created with the corresponding count matrices.

<details markdown="1">
<summary>Output files</summary>

- `summary/`
  - `*_(assignment|classification).csv`: Combined assignment/classification table from all used genetic- and hashing-based deconvolution tools.
  - `*_genetic.h5ad`: The RNA-seq count matrix with the assignments/classifications of all genetic tools saved in `.obs`.
  - `*_hashing.h5ad`: The hashing count matrix with the assignments/classifications of all hashing tools saved in `.obs`.
  - `*_genetic_and_hashing.h5mu`: Both `genetic.h5ad` and `hashing.h5ad` combined in a MuData object.

</details>

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
