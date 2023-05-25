# Output

This document describes the output produced by each process of the pipeline.

## Pipeline overview

<details><summary>Modes </summary>

- genetic: Genetics-based deconvolution workflow
- hashing: Hashing-based deconvolution workflow
- rescue: genetic + hashing + donor matching

</details>
<br>
<details><summary>Workflows </summary>
<details><summary>Hashing-based deconvolution (hash_demulti)</summary>

- Pre-processing
- Multiseq
- HTODemux
- HashedDrops
- DemuxEM
- HashSolo
- Solo
- Demuxmix
- GMM-Demux
- BFF

</details>

<details><summary>Genetics-based deconvolution (gene_demulti)</summary>

- Pre-processing: Samtools
- Variant-calling: freebayes
- Variant-filtering: BCFtools
- Variant-calling: cellsnp-lite
- Demuxlet
- Freemuxlet
- Vireo
- Souporcell
- scSplit

</details>
</details>

<br>

## Hashing-based deconvolution workflow

The output of hashing-based deconvolution workflow is saved in the folder `$projectDir/$params.outdir/$params.mode/hash_demulti`.

### Pre-processing

output directory: `preprocess/preprocess_[task_ID]`

- `${params.preprocessOut}.rds`: pre-processed data in an RDS object
- `params.csv`: specified parameters in the hashing pre-processing task

### HTODemux

output directory: `htodemux/htodemux_[task_ID]`

- `${params.assignmentOutHTO}_assignment_htodemux.csv`: the assignment of HTODemux
- `${params.assignmentOutHTO}_classification_htodemux.csv`: the classification of HTODemux as singlet, doublet and negative droplets
- `${params.objectOutHTO}.rds`: the result of HTODemux in an RDS object
- `params.csv`: specified parameters in the HTODemux task

Optionally:

- `ridge.jpeg`: a ridge plot showing the enrichment of selected HTOs
- `featureScatter.jpeg`: a scatter plot showing the signal of two selected HTOs
- `violinPlot.jpeg`: a violin plot showing selected features
- `tSNE.jpeg`: a 2D plot based on tSNE embedding of HTOs
- `heatMap.jpeg`: a heatmap of hashtag oligo signals across singlets, doublets and negative cells
- `visual_params.csv`: specified parameters for visualisation of the HTODemux result

### Multiseq

output directory: `multiseq/multiseq_[task_ID]`

- `${params.assignmentOutMulti}_res.csv`: the assignment of Multiseq
- `${params.objectOutMulti}.rds`: the result of Multiseq in an RDS object
- `params.csv`: specified parameters in the Multiseq task

### Demuxem

output directory: `demuxem/demuxem_[task_ID]`

- `${params.objectOutDemuxem}_demux.zarr.zip`: RNA expression matrix with demultiplexed sample identities in Zarr format
- `${params.objectOutDemuxem}.out.demuxEM.zarr.zip`: DemuxEM-calculated results in Zarr format, containing two datasets, one for HTO and one for RNA
- `${params.objectOutDemuxem}.ambient_hashtag.hist.pdf`: A histogram plot depicting hashtag distributions of empty droplets and non-empty droplets
- `${params.objectOutDemuxem}.background_probabilities.bar.pdf}`: A bar plot visualizing the estimated hashtag background probability distribution
- `${params.objectOutDemuxem}.real_content.hist.pdf`: A histogram plot depicting hashtag distributions of not-real-cells and real-cells as defined by total number of expressed genes in the RNA assay
- `${params.objectOutDemuxem}.rna_demux.hist.pdf`: This figure consists of two plots. The first one is a horizontal bar plot depicting the percentage of RNA barcodes with at least one HTO count. The second plot is a histogram plot depicting RNA UMI distribution for singlets, doublets and unknown cells.
- `${params..objectOutDemuxem}.gene_name.violin.pdf`: Violin plots depicting gender-specific gene expression across samples.
- `${params.objectOutDemuxem}_summary.csv`: the classification of Demuxem
- `${params.objectOutDemuxem}_obs.csv`: the assignment of Demuxem
- `params.csv`: specified parameters in the Demuxem task

Optionally:

- `{params.objectOutDemuxem}.{gene_name}.violin.pdf`: violin plots using specified gender-specific gene

### Solo

output directory: `solo/solo_[task_ID]`

- `${params.assignmentOutSolo}_res.csv`: the assignment of Solo
- `params.csv`: specified parameters in the Solo task

### HashSolo

output directory: `hashsolo/hashsolo_[task_ID]`

- `${params.assignmentOutHashSolo}_res.csv`: the assignment of HashSolo
- `${params.plotOutHashSolo}.jpg`: plot of HashSolo demultiplexing results for QC checks
- `params.csv`: specified parameters in the HashSolo task

### HashedDrops

output directory: `hashedDrops/hashedDrops_[task_ID]`

- `${params.objectOutEmptyDrops}.rds`: the result of emptyDrops in an RDS object
- `${params.assignmentOutEmptyDrops}.csv`: the result of emptyDrops in a csv file
- `plot_emptyDrops.png`: a diagnostic plot comparing the total count against the negative log-probability
- `${params.objectOutHashedDrops}.rds`: the result of hashedyDrops in an RDS object
- `${params.assignmentOutHashedDrops}_res.csv`: the assignment of HashSolo
- `${params.objectOutHashedDrops}_LogFC.png`: a diagnostic plot comparing the log-fold change between the second HTO's abundance and the ambient contamination
- `params.csv`: specified parameters in the HashedDrops task

### Demuxmix
output directory: `demuxmix/demuxmix_[task_ID]`

- `${params.assignmentOutDemuxmix}_assignment_demuxmix.csv`: the assignment and classification results produced by Demuxmix
- `params.csv`: specified parameters in the Demuxmix task

### GMM-Demux
output directory: `gmm_demux/gmm_demux_[task_ID]`

- `features.tsv.gz`: default content in the output folder are the non-MSM droplets (SSDs), stored in MTX format.
- `barcodes.tsv.gz`: default content in the output folder are the non-MSM droplets (SSDs), stored in MTX format.
- `matrix.mtx.gz`: default content in the output folder are the non-MSM droplets (SSDs), stored in MTX format.
- `GMM_full.csv`: The classification file containing the label of each droplet as well as the probability of the classification.
- `GMM_full.config`: Used to assign each classification to a donor using the numbers listed in the config file
- `gmm_demux_${task.index}_report.txt`: Specify the file to store summary report, produced only if GMM can find a viable solution that satisfies the droplet formation model
- `params.csv`: specified parameters in the GMM-Demux task

### BFF
output directory: `bff/bff_[task_ID]`

- `${params.assignmentOutBff}_assignment_demuxmix.csv`: the assignment and classification results produced by BFF
- `params.csv`: specified parameters in the BFF task

## Genetics-based deconvolution workflow

The output of genetics-based deconvolution workflow is saved in the folder `$projectDir/$params.outdir/$params.mode/gene_demulti`.

### Samtools

output directory: `samtools/samtools_[task_ID]`

- `filtered.bam`: processed BAM in a way that reads with any of following patterns be removed: read quality lower than 10, being unmapped segment, being secondary alignment, not passing filters, being PCR or optical duplicate, or being supplementary alignment
- `filtered.bam.bai`: index of filtered bam
- `no_dup.bam`: processed BAM after removing duplicated reads based on UMI
- `sorted.bam`: sorted BAM
- `sorted.bam.bai`: index of sorted BAM

### cellSNP-lite

- `cellSNP.base.vcf.gz`: a VCF file listing genotyped SNPs and aggregated AD & DP infomation (without GT)
- `cellSNP.samples.tsv`: a TSV file listing cell barcodes or sample IDs
- `cellSNP.tag.AD.mtx`: a file in mtx format, containing the allele depths of the alternative (ALT) alleles
- `cellSNP.tag.DP.mtx`: a file in mtx format, containing the sum of allele depths of the reference and alternative alleles (REF + ALT)
- `cellSNP.tag.OTH.mtx`: a file in mtx format, containing the sum of allele depths of all the alleles other than REF and ALT.
- `cellSNP.cells.vcf.gz`: a VCF file listing genotyped SNPs and AD & DP & genotype (GT) information for each cell or sample
- `params.csv`: specified parameters in the cellsnp-lite task

### Freebayes

- `${region}_${vcf_freebayes}`: a VCF file containing variants called from mixed samples in the given chromosome region

### Bcftools

output directory: `bcftools/bcftools_[task_ID]`

- `total_chroms.vcf`: a VCF containing variants from all chromosomes
- `sorted_total_chroms.vcf`: sorted VCF file
- `filtered_sorted_total_chroms.vcf`: sorted VCF file containing variants with a quality score > 30

### Demuxlet

output directory: `demuxlet/demuxlet_[task_ID]`

- `{demuxlet_out}.best`: result of demuxlet containing the best guess of the sample identity, with detailed statistics to reach to the best guess
- `params.csv`: specified parameters in the Demuxlet task

Optionally:

- `{demuxlet_out}.cel`: contains the relation between numerated barcode ID and barcode. Also, it contains the number of SNP and number of UMI for each barcoded droplet.
- `{demuxlet_out}.plp`: contains the overlapping SNP and the corresponding read and base quality for each barcode ID.
- `{demuxlet_out}.umi`: contains the position covered by each umi
- `{demuxlet_out}.var`: contains the position, reference allele and allele frequency for each SNP.

### Freemuxlet

output directory: `freemuxlet/freemuxlet_[task_ID]`

- `{freemuxlet_out}.clust1.samples.gz`: contains the best guess of the sample identity, with detailed statistics to reach to the best guess.
- `{freemuxlet_out}.clust1.vcf.gz`: VCF file for each sample inferred and clustered from freemuxlet
- `{freemuxlet_out}.lmix`: contains basic statistics for each barcode
- `params.csv`: specified parameters in the Freemuxlet task

Optionally:

- `{freemuxlet_out}.cel`: contains the relation between numerated barcode ID and barcode. Also, it contains the number of SNP and number of UMI for each barcoded droplet.
- `{freemuxlet_out}.plp`: contains the overlapping SNP and the corresponding read and base quality for each barcode ID.
- `{freemuxlet_out}.umi`: contains the position covered by each umi
- `{freemuxlet_out}.var`: contains the position, reference allele and allele frequency for each SNP.
- `{freemuxlet_out}.clust0.samples.gz`: contains the best sample identity assuming all droplets are singlets
- `{freemuxlet_out}.clust0.vcf.gz}`: VCF file for each sample inferred and clustered from freemuxlet assuming all droplets are singlets
- `{freemuxlet_out}.ldist.gz`: contains the pairwise Bayes factor for each possible pair of droplets

### Vireo

output directory: `vireo/vireo_[task_ID]`

- `donor_ids.tsv`: assignment of Vireo with detailed statistics
- `summary.tsv`: summary of assignment
- `prob_singlet.tsv.gz`: contains probability of classifing singlets
- `prob_doublet.tsv.gz`: contains probability of classifing doublets
- `GT_donors.vireo.vcf.gz`: contains estimated donor genotypes
- `filtered_variants.tsv`: a minimal set of discriminatory variants
- `GT_barcodes.png`: a figure for the identified genotype barcodes
- `fig_GT_distance_estimated.pdf`: a plog showing estimated genotype distance
- `_log.txt`: vireo log file
- `params.csv`: specified parameters in the Vireo task

### scSplit

output directory: `scSplit/scsplit_[task_ID]`

- `alt_filtered.csv`: count matrix of alternative alleles
- `ref_filtered.csv`: count matrix of reference alleles
- `scSplit_result.csv`: barcodes assigned to each of the N+1 cluster (N singlets and 1 doublet cluster), doublet marked as DBL-<n> (n stands for the cluster number), e.g SNG-0 means the cluster 0 is a singlet cluster.
- `scSplit_dist_matrix.csv`: the ALT allele Presence/Absence (P/A) matrix on distinguishing variants for all samples as a reference in assigning sample to clusters, NOT including the doublet cluster, whose sequence number would be different every run (please pay enough attention to this)
- `scSplit_dist_variants.txt`: the distinguishing variants that can be used to genotype and assign sample to clusters
- `scSplit_PA_matrix.csv`: the full ALT allele Presence/Absence (P/A) matrix for all samples, NOT including the doublet cluster, whose sequence number would be different every run (please pay enough attention to this)
- `scSplit_P_s_c.csv`: the probability of each cell belonging to each sample
- `scSplit.log`: log file containing information for current run, iterations, and final Maximum Likelihood and doublet sample
- `params.csv`: specified parameters in the scSplit task

### Souporcell

output directory: `souporcell/souporcell_[task_ID]`

- `alt.mtx`: count matrix of alternative alleles
- `ref.mtx`: count matrix of reference alleles
- `clusters.tsv`: assignment of Souporcell with the cell barcode, singlet/doublet status, cluster, log_loss_singleton, log_loss_doublet, followed by log loss for each cluster.
- `cluster_genotypes.vcf`: VCF with genotypes for each cluster for each variant in the input vcf from freebayes
- `ambient_rna.txt`: contains the ambient RNA percentage detected
- `params.csv`: specified parameters in the Souporcell task

## Merging results

After each demultiplexing workflow, the pipeline will generate some TSV files to summarize the results in the folder `$projectDir/$params.outdir/$params.mode/[workflow]/[workflow]_summary`.

- `[method]_classification.csv`: classification of all trials for a given method
  | Barcode | multiseq_1 | multiseq_2 | ... |
  |:---------: |:----------: |:----------: |:---: |
  | barcode-1 | singlet | singlet | ... |
  | barcode-2 | doublet | negative | ... |
  | ... | ... | ... | ... |
- `[method]_assignment.csv`: assignment of all trials for a given method
  | Barcode | multiseq_1 | multiseq_2 | ... |
  |:---------: |:----------: |:----------: |:---: |
  | barcode-1 | donor-1 | donor-2 | ... |
  | barcode-2 | doublet | negative | ... |
  | ... | ... | ... | ... |
- `[method]_params.csv`: specified paramters of all trials for a given method
  | Argument | Value |
  | :---------: | :----------: |
  | seuratObejctPath | Path |
  | quantile | 0.7 |
  | ... | ... |
- `[workflow]_classification_all.csv`: classification of all trials across different methods
  | Barcode | multiseq_1 | htodemux_1 | ... |
  |:---------: |:----------: |:----------: |:---: |
  | ... | ... | ... | ... |
- `[workflow]_assignment_all.csv`: save the assignment of all trials across different methods
  | Barcode | multiseq_1 | htodemux_1 | ... |
  |:---------: |:----------: |:----------: |:---: |
  | ... | ... | ... | ... |

- In the `rescue` mode, the pipeline merges the results of hashing and genetic demultiplexing tools into and `assignment_all_genetic_and_hash.csv` in the `$projectDir/$params.outdir/$params.mode/summary` folder.

## Donor matching

Output directory: `$projectDir/$params.outdir/$params.mode/donor_match`

- folder`[method1]_[task_ID]_vs_[method2]_[task_ID]` with:
  - `correlation_res.csv`: correlation scores of donor matching
  - `concordance_heatmap.png`: a heatmap visualising the the correlation scores
  - `donor_match.csv`: a map between hashtag and donor identity.
- For the optimal match `best_method1` and `best_method2` among all trials, the pipeline generates new assignment file by mapping hashtags of `best_method2` to `best_method1` in folder `donor_match/donor_match_[best_method1]_[best_method2]`:
  - `all_assignment_after_match.csv`: assignment of all cell barcodes after donor matching
  - `intersect_assignment_after_match.csv`: assignment of joint singlets after donor matching
  - Optionally, if `best_method1` is `vireo` and identification of donor-specific or discriminatory variants is enabled:
    - `donor_specific_variants.csv`: a list of donor-specific variants
    - `donor_specific_variants_upset.png`: An upset plot showing the number of donor-specific variants
    - `donor_genotype_subset_by_default_matched.vcf`: Donor genotypes of donor-specific variants
    - `donor_genotype_subset_by_vireo.vcf`: Donor genotypes of a set of discriminatory variants filtered by Vireo
