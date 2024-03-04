# Running hadge on multiple samples

The pipeline is able to run on multiple samples.
In this scenario, the shared parameters for input data are retrieved from a sample sheet using `params.multi_sample`, which is set to `null` by default.

## **Sample sheet**

- The sample sheet should contain a column called `sampleId` for unique sample IDs assigned to each sample.
- The sample sheet (example file see the Resources section below) must contain different columns depending on the mode and methods you want to run.

  - hashing mode:
    | sampleId | rna_matrix_raw | rna_matrix_filtered | hto_matrix_raw | hto_matrix_filtered |
    |----------|----------------|---------------------|----------------|---------------------|
    | sample1 | | | | |
    | sample2 | | | | |
  - genetic mode: Set the value to "None" if the input data, for example, vcf_donor, is not available, similar to the single-sample mode. Do not forget to include the columns for HTO and RNA count matrices if `params.generate_anndata` or `params.generate_mudata` is enabled.
    | sampleId | bam | bam_index | barcodes | nsample | celldata | vcf_mixed | vcf_donor |
    |----------|-----|-----------|----------|---------|----------|-----------|-----------|
    | sample1 | | | | | | | |
    | sample2 | | | | | | | |
  - rescue mode:
    | sampleId | rna_matrix_raw | rna_matrix_filtered | hto_matrix_raw | hto_matrix_filtered | bam | bam_index | barcodes | nsample | celldata | vcf_mixed | vcf_donor |
    |----------|----------------|---------------------|----------------|---------------------|-----|-----------|----------|---------|----------|-----------|-----------|
    | sample1 | | | | | | | | | | | |
    | sample2 | | | | | | | | | | | |

- The remaining parameters for each process are specified in the nextflow.config file, just like when demultiplexing a single sample.
- There is a distinction between running on a single sample and running on multiple samples. When processing multiple samples, the pipeline only permits a single value for each process parameter, whereas in the case of a single sample, multiple values separated by commas are allowed.

## **Output**

When running the pipeline on multiple samples, the pipeline output will be found in the folder `"$projectDir/$params.outdir/$sampleId/$params.mode`.

## **Resources**

There is an [example sample sheet](../../test_data/multi_sample_input.csv) for `multi_sample` mode.
