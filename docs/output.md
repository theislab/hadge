# nf-core/hadge: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

### Genetic-based deconvolution

<details markdown="1">
<summary>Output files</summary>

- `genetic/`
  - `popscle/demuxlet/`: See [modules/popscle_demuxlet](https://nf-co.re/modules/popscle_demuxlet/) for detailed output information.
  - `popscle/freemuxlet/`: See [modules/popscle_freemuxlet](https://nf-co.re/modules/popscle_freemuxlet/) for detailed output information.
  - `souporcell/`: See [modules/souporcell](https://nf-co.re/modules/souporcell/) for detailed output information.
  - `vireo/`: See [modules/vireo](https://nf-co.re/modules/vireo/) for detailed output information.
- `genetic/summary/`
  - `/*/*_genetic_summary_(assignment|classification).csv`: Summary of all assigned/classified cells from each genetic-based deconvolution tool, merged into a single table.
  - `/*/*_genetic_overview_(assignment|classification).csv`: This table summarizes each genetic-based deconvolution tool (before merging) by reporting its total barcode count, the number of barcodes it shares with every other method, and the counts of each donor label or classification category (e.g., `0`, `1`, `singlet`, `doublet`, `negative`).

</details>

This subworkflow runs when the `--mode` is set to `genetic` or `rescue`. It saves the results of each genetic-based deconvolution specified with `--genetic_tools` in a folder named after the corresponding tool. Additionally, this step summarizes all deconvolution results in a single summary table. The assignments CSV provides, for each cell barcode, one of the following categories:

| Category      | Description                                    |
| ------------- | ---------------------------------------------- |
| Cluster label | A cluster label (e.g., `0`, `1`, …)            |
| `doublet`     | More than one cluster is assigned to this cell |
| `negative`    | All other cases (e.g., undetected cells)       |

Classification works similarly, except that cluster labels are replaced by the label `singlet`.

### Hashing-based deconvolution

<details markdown="1">
<summary>Output files</summary>

- `hashing/`
  - `bff/`: See [modules/bff](https://nf-co.re/modules/bff/) for detailed output information.
  - `demuxem/`: See [modules/demuxem](https://nf-co.re/modules/bff/) for detailed output information.
  - `gmm-demux/`: See [modules/gmmdemux](https://nf-co.re/modules/gmmdemux/) for detailed output information.
  - `hasheddrops/`: See [modules/hasheddrops](https://nf-co.re/modules/hasheddrops/) for detailed output information.
  - `hashsolo/`: See [modules/scanpy_hashsolo](https://nf-co.re/modules/scanpy_hashsolo/) for detailed output information.
  - `htodemux/`: See [modules/htodemux](https://nf-co.re/modules/htodemux/) for detailed output information.
  - `htodemux/visualization/`: Visualizations of htodemux results.
  - `multiseq/`: See [modules/multiseqdemux](https://nf-co.re/modules/multiseqdemux/) for detailed output information.
- `hashing/summary/`
  - `/*/*_hashing_summary_(assignment|classification).csv`: Summary of all assigned/classified cells from each hashing-based deconvolution tool, merged into a single table.
  - `/*/*_hashing_overview_(assignment|classification).csv`: This table summarizes each hashing-based deconvolution tool (before merging) by reporting its total barcode count, the number of barcodes it shares with every other method, and the counts of each donor label or classification category (e.g., `HTO-1`, `HTO-2`, `singlet`, `doublet`, `negative`).

</details>

This subworkflow runs when the `--mode` is set to `hashing` or `rescue`. It saves the results of each hashing-based deconvolution specified with `--hash_tools` in a folder named after the corresponding tool. Additionally, this step summarizes all deconvolution results in a single summary table. The assignments CSV provides, for each cell barcode, one of the following categories:

| Category    | Description                                                       |
| ----------- | ----------------------------------------------------------------- |
| Donor label | A HTO label that identifies the donor (e.g., `HTO-1`, `HTO-2`, …) |
| `doublet`   | More than one donor is assigned to this cell                      |
| `negative`  | All other cases (e.g., undetected cells)                          |

Classification works similarly, except that cluster labels are replaced by the label `singlet`.

### Donor matching

<details markdown="1">
<summary>Output files</summary>

- `donor_match/`
  - `*_best_all_assignment_after_match.csv`: assignment of all cell barcodes based on the donor matching of the optimal match
  - `*_best_donor_match.csv`: a map between hashtags and donor identities based on the donor matching of the optimal match
  - `*_best_intersect_assignment_after_match.csv`: assignment of joint singlets based on the donor matching of the optimal match
  - `*_score_record.csv`: a CSV file storing the matching score and the number of matched donors for each method pair
- `donor_match/[method1]_vs_[method2]/`
  - `*_all_assignment_after_match.csv`: assignment of all cell barcodes after donor matching
  - `*_concordance_heatmap.png`: a heatmap visualising the the correlation scores
  - `*_correlation_res.csv`: correlation scores of donor matching
  - `*_donor_match.csv`: a map between hashtag and donor identity.
  - `*_intersect_assignment_after_match.csv`: assignment of joint singlets after donor matching

</details>

For each hashing–genetic deconvolution method pair, pairwise Pearson correlations are computed between binarized cell assignment vectors to match donors.
A matching score is obtained by summing the non-negative correlations and dividing by the number of expected donors, providing a measure of agreement between methods.

Donor matching is executed only when `--match_donor` is enabled (default).
Full output is produced in `rescue` and `donor_match` mode, while other modes generate a reduced set of files, as matching is restricted to comparisons within the same method type (genetic-to-genetic and hashing-to-hashing) and not across types.

### Find variants

<details markdown="1">
<summary>Output files</summary>

- `find_variants/`
  - `*_all_representative_variants.csv`: a list of representative variants from all donors
  - `*_donor_specific_variants.csv`: a list of donor-specific variants
  - `*_donor_specific_variants_upset.png`: an upset plot showing the number of donor-specific variants
  - `*_vireo_variants.csv`: a list of discriminatory variants filtered by Vireo
- `donor_match/[hto_name]/`
  - `*_informative_variants.csv`: informative variants for the donor indicated in the file name
  - `*_matched_gt.csv`: genotype of the donor indicated in the file name
  - `*_unmatched_gt.csv`: genotype of all other donors
- `donor_match/subset_gt_donors/`
  - `_donor_specific.vcf.gz`: Donor genotypes of donor-specific variants
  - `_vireo.vcf.gz`: Donor genotypes of a set of discriminatory variants filtered by Vireo

</details>

Find variants runs only when `--match_donor` is enabled (default) and the workflow is in `rescue` or `donor_match` mode.
It outputs two sets of variants: informative variants from vireo and donor specific variants as described in the [paper](https://link.springer.com/article/10.1186/s13059-024-03249-z#:~:text=Vireo%20is%20implemented,others%20during%20deconvolution.).
A subset VCF of donor genotypes (produced by vireo) is then generated.

### Summary

<details markdown="1">
<summary>Output files</summary>

- `summary/`
  - `*_(assignment|classification).csv`: Combined assignment/classification table from all used genetic- and hashing-based deconvolution tools.
  - `*_genetic.h5ad`: The RNA-seq count matrix with the assignments/classifications of all genetic tools saved in `.obs`. If the `barcodes.tsv` from the samplesheet input contains barcodes not present in the RNA-seq count matrix, those barcodes will be dropped when joining the results from the genetic tools. Therefore, it is recommended that both files share the same set of barcodes.
  - `*_hashing.h5ad`: The hashing count matrix with the assignments/classifications of all hashing tools saved in `.obs`.
  - `*_genetic_and_hashing.h5mu`: Both `genetic.h5ad` and `hashing.h5ad` combined in a MuData object.

</details>

When running the pipeline in `rescue` mode, the combined assignment/classification table from all used genetic- and hashing-based deconvolution tools will be generated here.
Additionally, an AnnData/MuData object is created with the corresponding count matrices.

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
