#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process freebayes{
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti/freebayes", mode: 'copy'

    input:
        file bam_freebayes
        file bai_freebayes
        val bam_list_freebayes
        val stdin_freebayes
        val ref_freebayes
        val ref_index_freebayes
        val targets_freebayes
        each region_freebayes
        val samples_freebayes
        val populations_freebayes
        val cnv_map_freebayes
        val vcf_freebayes
        val gvcf_freebayes
        val gvcf_chunk_freebayes
        val gvcf_dont_use_chunk_freebayes
        val variant_input_freebayes
        val only_use_input_alleles_freebayes
        val haplotype_basis_alleles_freebayes
        val report_all_haplotype_alleles_freebayes
        val report_monomorphic_freebayes
        val pvar_freebayes
        val strict_vcf_freebayes
        val theta_freebayes
        val ploidy_freebayes
        val pooled_discrete_freebayes
        val pooled_continuous_freebayes
        val use_reference_allele_freebayes
        val reference_quality_freebayes
        val no_snps_freebayes
        val no_indels_freebayes
        val no_mnps_freebayes
        val no_complex_freebayes
        val use_best_n_alleles_freebayes
        val haplotype_length_freebayes
        val min_repeat_size_freebayes
        val min_repeat_entropy_freebayes
        val no_partial_observations_freebayes
        val dont_left_align_indels_freebayes
        val use_duplicate_reads_freebayes
        val min_mapping_quality_freebayes
        val min_base_quality_freebayes
        val min_supporting_allele_qsum_freebayes
        val min_supporting_mapping_qsum_freebayes
        val mismatch_base_quality_threshold_freebayes
        val read_mismatch_limit_freebayes
        val read_max_mismatch_fraction_freebayes
        val read_snp_limit_freebayes
        val read_indel_limit_freebayes
        val standard_filters_freebayes
        val min_alternate_fraction_freebayes
        val min_alternate_count_freebayes
        val min_alternate_qsum_freebayes
        val min_alternate_total_freebayes
        val min_coverage_freebayes
        val max_coverage_freebayes
        val no_population_priors_freebayes
        val hwe_priors_off_freebayes
        val binomial_obs_priors_off_freebayes
        val allele_balance_priors_off_freebayes
        val observation_bias_freebayes
        val base_quality_cap_freebayes
        val prob_contamination_freebayes
        val legacy_gls_freebayes
        val contamination_estimates_freebayes
        val report_genotype_likelihood_max_freebayes
        val genotyping_max_iterations_freebayes
        val genotyping_max_banddepth_freebayes
        val posterior_integration_limits_freebayes
        val exclude_unobserved_genotypes_freebayes
        val genotype_variant_threshold_freebayes
        val use_mapping_quality_freebayes
        val harmonic_indel_quality_freebayes
        val read_dependence_factor_freebayes
        val genotype_qualities_freebayes
        val debug_freebayes
        val dd_freebayes
        

    output:
        path "*.vcf"

    script:
    def bam_list = bam_list_freebayes != 'no_bam_list' ? "--bam-list ${bam_list}":''
    def stdin = stdin_freebayes != 'False' ? "--stdin" : ''
    def targets = targets_freebayes != 'no_targets' ? "--targets ${targets_freebayes}" : ''
    def region = region_freebayes != 'False' ? "--region ${region_freebayes}" : ''
    def samples = samples_freebayes != 'no_samples' ? "--samples ${samples_freebayes}" : ''
    def populations = populations_freebayes != 'no_populations' ? "--populations ${populations}":''
    def cnv_map = cnv_map_freebayes != 'no_cnv_map' ? "--cnv-map ${cnv_map_freebayes}" : ''
 
    def gvcf = gvcf_freebayes !='False' ? "--gvcf":''
    def gvcf_chunk = gvcf_chunk_freebayes !='False' ? "--gvcf-chunk ${gvcf_chunk_freebayes}":''
    def gvcf_dont_use_chunk = gvcf_dont_use_chunk_freebayes != 'False' ?"--gvcf-dont-use-chunk ${gvcf_dont_use_chunk_freebayes}":''
    def variant_input = variant_input_freebayes != 'no_variant_input' ? "--variant-input ${variant_input_freebayes}" : ''
    def only_use_input_alleles = only_use_input_alleles_freebayes != 'False' ? "--only-use-input-alleles" : ''
    def haplotype_basis_alleles = haplotype_basis_alleles_freebayes != 'no_haplotype_basis_alleles' ? "--haplotype-basis-alleles ${haplotype_basis_alleles_freebayes}":''
    def report_all_haplotype_alleles = report_all_haplotype_alleles_freebayes !='False' ? "--report-all-haplotype-alleles":''
    def report_monomorphic = report_monomorphic_freebayes !='False' ? "--report-monomorphic":''
    def pvar = "--pvar ${pvar_freebayes}"
    def strict_vcf = strict_vcf_freebayes!='False' ? "--strict-vcf":''
    def theta = "--theta ${theta_freebayes}"
    def ploidy = "--ploidy ${ploidy_freebayes}"
    def pooled_discrete = pooled_discrete_freebayes!='False' ? "--pooled-discrete":''
    def pooled_continuous = pooled_continuous_freebayes !='False' ? "--pooled-continuous":''

    def use_reference_allele = use_reference_allele_freebayes != 'False' ? "--use-reference-allele" : ''
    def reference_quality = "--reference-quality ${reference_quality_freebayes}"

    def no_snps = no_snps_freebayes !='False' ? "--no-snps" :''
    def no_indels = no_indels_freebayes !='False' ? "--no-indels" :''
    def no_mnps = no_mnps_freebayes !='False' ? "--no-mnps" :''
    def no_complex = no_complex_freebayes !='False' ? "--no-complex" :''

    def use_best_n_alleles = "--use-best-n-alleles ${use_best_n_alleles_freebayes}"
    def haplotype_length = haplotype_length_freebayes = "--haplotype-length ${haplotype_length_freebayes}"
    def min_repeat_size = "--min-repeat-size ${min_repeat_size_freebayes}"
    def min_repeat_entropy = "--min-repeat-entropy ${min_repeat_entropy_freebayes}"
    def no_partial_observations = no_partial_observations_freebayes !='False' ? "--no-partial-observations":''

    def dont_left_align_indels = dont_left_align_indels_freebayes != 'False' ? "--dont-left-align-indels" : ''

    def use_duplicate_reads = use_duplicate_reads_freebayes != 'False' ? "--use-duplicate-reads" : ''
    def min_mapping_quality = "--min-mapping-quality ${min_mapping_quality_freebayes}"
    def min_base_quality = "--min-base-quality ${min_base_quality_freebayes}"
    def min_supporting_allele_qsum = "--min-supporting-allele-qsum ${min_supporting_allele_qsum_freebayes}"
    def min_supporting_mapping_qsum = "--min-supporting-mapping-qsum ${min_supporting_mapping_qsum_freebayes}"
    def mismatch_base_quality_threshold = "--mismatch-base-quality-threshold ${mismatch_base_quality_threshold_freebayes}"
    def read_mismatch_limit = read_mismatch_limit_freebayes != 'False' ? "-read-mismatch-limit ${read_mismatch_limit_freebayes}":''
    def read_max_mismatch_fraction = "--read-max-mismatch-fraction ${read_max_mismatch_fraction_freebayes}"
    def read_snp_limit = read_snp_limit_freebayes != 'False' ? "--read-snp-limit ${read_snp_limit_freebayes}":''
    def read_indel_limit = read_indel_limit_freebayes != 'False' ? "--read-indel-limit ${read_indel_limit_freebayes}" :''
    def standard_filters = standard_filters_freebayes!='False' ? "--standard-filters":''
    def min_alternate_fraction = "--min-alternate-fraction ${min_alternate_fraction_freebayes}"
    def min_alternate_count = "--min-alternate-count ${min_alternate_count_freebayes}"
    def min_alternate_qsum = "--min-alternate-qsum ${min_alternate_qsum_freebayes}"
    def min_alternate_total = "--min-alternate-total ${min_alternate_total_freebayes}"
    def min_coverage = "--min-coverage ${min_coverage_freebayes}"
    def max_coverage = max_coverage_freebayes !='False' ? "--max-coverage ${max_coverage_freebayes}":''

    def no_population_priors = no_population_priors_freebayes != 'False' ? "--no-population-priors" : ''
    def hwe_priors_off = hwe_priors_off_freebayes != 'False' ? "--hwe-priors-off" : ''
    def binomial_obs_priors_off = binomial_obs_priors_off_freebayes != 'False' ? "--binomial-obs-priors-off" : ''
    def allele_balance_priors_off = allele_balance_priors_off_freebayes != 'False' ? "--allele-balance-priors-off" : ''
    def observation_bias = observation_bias_freebayes!= 'no_observation_bias' ? "--observation-bias ${observation_bias_freebayes}":''
    def base_quality_cap = base_quality_cap_freebayes !='False' ? "--base-quality-cap ${base_quality_cap_freebayes}":''
    def prob_contamination = "--prob-contamination ${prob_contamination_freebayes}"
    def legacy_gls = legacy_gls_freebayes !='False' ? "--legacy-gls" :''
    def contamination_estimates = contamination_estimates_freebayes != 'no_contamination_estimates' ? "--contamination-estimates ${contamination_estimates_freebayes}":''
    
    def report_genotype_likelihood_max = report_genotype_likelihood_max_freebayes !='False' ? "--report-genotype-likelihood-max":''
    def genotyping_max_iterations = "--genotyping-max-iterations ${genotyping_max_iterations_freebayes}"
    def genotyping_max_banddepth = "--genotyping-max-banddepth ${genotyping_max_banddepth_freebayes}"
    def posterior_integration_limits = "--posterior-integration-limits ${posterior_integration_limits_freebayes}"
    def exclude_unobserved_genotypes = exclude_unobserved_genotypes_freebayes != 'False' ? "--exclude-unobserved-genotypes" : ''
    def genotype_variant_threshold = genotype_variant_threshold_freebayes != 'False' ? "--genotype-variant-threshold ${genotype_variant_threshold_freebayes}":''
    def use_mapping_quality = use_mapping_quality_freebayes != 'False' ? "--use-mapping-quality" : ''
    def harmonic_indel_quality = harmonic_indel_quality_freebayes !='False' ? "--harmonic-indel-quality" :''
    def read_dependence_factor = "--read-dependence-factor ${read_dependence_factor_freebayes}"
    def genotype_qualities = genotype_qualities_freebayes !='False' ? "--genotype-qualities" :''

    def debug = debug_freebayes != 'False' ? "--debug" : ''
    def dd = dd_freebayes != 'False' ? "-dd" : ''

    """
    freebayes ${bam_freebayes} ${bam_list} $stdin -f ${ref_freebayes} $targets $region $samples $populations ${cnv_map} \
    -v ${region_freebayes}_${vcf_freebayes} $gvcf ${gvcf_chunk} ${gvcf_dont_use_chunk} ${variant_input} ${only_use_input_alleles} ${haplotype_basis_alleles} ${report_all_haplotype_alleles} ${report_monomorphic} $pvar ${strict_vcf} \
    $theta $ploidy ${pooled_discrete} ${pooled_continuous} \
    ${use_reference_allele} ${reference_quality} ${no_snps} ${no_indels} ${no_mnps} ${no_complex} ${use_best_n_alleles} ${haplotype_length} ${min_repeat_size} ${min_repeat_entropy} ${no_partial_observations} ${dont_left_align_indels} \
    ${use_duplicate_reads} ${min_mapping_quality} ${min_base_quality} ${min_supporting_allele_qsum} ${min_supporting_mapping_qsum} ${mismatch_base_quality_threshold} \
    ${read_mismatch_limit} ${read_max_mismatch_fraction} ${read_snp_limit} ${read_indel_limit} ${standard_filters} ${min_alternate_fraction} ${min_alternate_count} ${min_alternate_qsum} ${min_alternate_total} ${min_coverage} ${max_coverage} \
    ${no_population_priors} ${hwe_priors_off} ${binomial_obs_priors_off} ${allele_balance_priors_off} \
    ${observation_bias} ${base_quality_cap} ${prob_contamination} ${legacy_gls} ${contamination_estimates} \
    ${report_genotype_likelihood_max} ${genotyping_max_iterations} ${genotyping_max_banddepth} ${posterior_integration_limits} ${exclude_unobserved_genotypes} ${genotype_variant_threshold} \
    ${use_mapping_quality} ${harmonic_indel_quality} ${read_dependence_factor} ${genotype_qualities} $debug $dd
  
    """
}


workflow variant_freebayes{
    take:
        bam_freebayes
        bai_freebayes
        region_freebayes
    main:
        bam_list_freebayes = channel.value(params.bam_list)
        stdin_freebayes = Channel.value(params.stdin)
        fasta_reference = Channel.value(params.fasta)
        fasta_reference_index = Channel.value(params.fasta_index)
        targets_freebayes = Channel.value(params.targets)
        samples_freebayes = Channel.value(params.samples)
        populations_freebayes = channel.value(params.populations)
        cnv_map_freebayes = Channel.value(params.cnv_map)
        vcf_freebayes = channel.value(params.vcf_freebayes)
        gvcf_freebayes = channel.value(params.gvcf)
        gvcf_chunk_freebayes = channel.value(params.gvcf_chunk)
        gvcf_dont_use_chunk_freebayes = channel.value(params.gvcf_dont_use_chunk)
        variant_input_freebayes = Channel.value(params.variant_input)
        only_use_input_alleles_freebayes = Channel.value(params.only_use_input_alleles)
        haplotype_basis_alleles_freebayes = channel.value(params.haplotype_basis_alleles)
        report_all_haplotype_alleles_freebayes = channel.value(params.report_all_haplotype_alleles)
        report_monomorphic_freebayes = channel.value(params.report_monomorphic)
        pvar_freebayes = channel.value(params.pvar)
        strict_vcf_freebayes = channel.value(params.strict_vcf)
        
        theta_freebayes = channel.value(params.theta)
        ploidy_freebayes = channel.value(params.ploidy)
        pooled_discrete_freebayes = channel.value(params.pooled_discrete)
        pooled_continuous_freebayes = channel.value(params.pooled_continuous)
        use_reference_allele_freebayes = channel.value(params.use_reference_allele)
        reference_quality_freebayes = channel.value(params.reference_quality)
        no_snps = channel.value(params.no_snps)
        no_indels = channel.value(params.no_indels)
        no_mnps = channel.value(params.no_mnps)
        no_complex = channel.value(params.no_complex)
        use_best_n_alleles_freebayes = channel.value(params.use_best_n_alleles)
        haplotype_length_freebayes = channel.value(params.haplotype_length)
        min_repeat_size_freebayes = channel.value(params.min_repeat_size)
        min_repeat_entropy_freebayes = channel.value(params.min_repeat_entropy)
        no_partial_observations_freebayes = channel.value(params.no_partial_observations)
        
        dont_left_align_indels_freebayes = channel.value(params.dont_left_align_indels)
        use_duplicate_reads_freebayes = channel.value(params.use_duplicate_reads)
        min_mapping_quality_freebayes = channel.value(params.min_mapping_quality)
        min_base_quality_freebayes = channel.value(params.min_base_quality)
        min_supporting_allele_qsum_freebayes = channel.value(params.min_supporting_allele_qsum)
        min_supporting_mapping_qsum_freebayes = channel.value(params.min_supporting_mapping_qsum)
        mismatch_base_quality_threshold_freebayes = channel.value(params.mismatch_base_quality_threshold)
        read_mismatch_limit_freebayes = channel.value(params.read_mismatch_limit)
        read_max_mismatch_fraction_freebayes = channel.value(params.read_max_mismatch_fraction)
        read_snp_limit_freebayes = channel.value(params.read_snp_limit)
        read_indel_limit_freebayes = channel.value(params.read_indel_limit)
        standard_filters_freebayes = channel.value(params.standard_filters)
        min_alternate_fraction_freebayes = channel.value(params.min_alternate_fraction)
        min_alternate_count_freebayes = channel.value(params.min_alternate_count)
        min_alternate_qsum_freebayes = channel.value(params.min_alternate_qsum)
        min_alternate_total_freebayes = channel.value(params.min_alternate_total)
        min_coverage_freebayes = channel.value(params.min_coverage)
        max_coverage_freebayes = channel.value(params.max_coverage)
        no_population_priors_freebayes = channel.value(params.no_population_priors)
        
        hwe_priors_off_freebayes = channel.value(params.hwe_priors_off)
        binomial_obs_priors_off_freebayes = channel.value(params.binomial_obs_priors_off)
        allele_balance_priors_off_freebayes = channel.value(params.allele_balance_priors_off)
        observation_bias_freebayes = channel.value(params.observation_bias)
        base_quality_cap_freebayes = channel.value(params.base_quality_cap)
        prob_contamination_freebayes = channel.value(params.prob_contamination)
        legacy_gls_freebayes = channel.value(params.legacy_gls)
        contamination_estimates_freebayes = channel.value(params.contamination_estimates)
        
        report_genotype_likelihood_max_freebayes = channel.value(params.report_genotype_likelihood_max)
        genotyping_max_iterations_freebayes = channel.value(params.genotyping_max_iterations)
        genotyping_max_banddepth_freebayes = channel.value(params.genotyping_max_banddepth)
        posterior_integration_limits_freebayes = channel.value(params.posterior_integration_limits)
        exclude_unobserved_genotypes_freebayes = channel.value(params.exclude_unobserved_genotypes)
        genotype_variant_threshold_freebayes = channel.value(params.genotype_variant_threshold)
        use_mapping_quality_freebayes = channel.value(params.use_mapping_quality)
        harmonic_indel_quality_freebayes = channel.value(params.harmonic_indel_quality)
        read_dependence_factor_freebayes = channel.value(params.read_dependence_factor)
        genotype_qualities_freebayes = channel.value(params.genotype_qualities)
        debug_freebayes = channel.value(params.debug)
        dd_freebayes = channel.value(params.dd)
        

        freebayes(bam_freebayes, bai_freebayes, bam_list_freebayes, stdin_freebayes, fasta_reference, fasta_reference_index, targets_freebayes, region_freebayes, samples_freebayes, populations_freebayes, cnv_map_freebayes, \
        vcf_freebayes, gvcf_freebayes, gvcf_chunk_freebayes, gvcf_dont_use_chunk_freebayes, variant_input_freebayes, only_use_input_alleles_freebayes, haplotype_basis_alleles_freebayes, report_all_haplotype_alleles_freebayes, report_monomorphic_freebayes, pvar_freebayes, strict_vcf_freebayes, \
        theta_freebayes, ploidy_freebayes, pooled_discrete_freebayes, pooled_continuous_freebayes, use_reference_allele_freebayes, reference_quality_freebayes, \
        no_snps, no_indels, no_mnps, no_complex, use_best_n_alleles_freebayes, haplotype_length_freebayes, min_repeat_size_freebayes, min_repeat_entropy_freebayes, no_partial_observations_freebayes, \
        dont_left_align_indels_freebayes, use_duplicate_reads_freebayes, min_mapping_quality_freebayes, min_base_quality_freebayes, min_supporting_allele_qsum_freebayes,\
        min_supporting_mapping_qsum_freebayes, mismatch_base_quality_threshold_freebayes, read_mismatch_limit_freebayes, read_max_mismatch_fraction_freebayes,
        read_snp_limit_freebayes, read_indel_limit_freebayes, standard_filters_freebayes, min_alternate_fraction_freebayes, min_alternate_count_freebayes, min_alternate_qsum_freebayes, min_alternate_total_freebayes, min_coverage_freebayes, max_coverage_freebayes, \
        no_population_priors_freebayes, hwe_priors_off_freebayes, binomial_obs_priors_off_freebayes, allele_balance_priors_off_freebayes, observation_bias_freebayes, base_quality_cap_freebayes, prob_contamination_freebayes, legacy_gls_freebayes, contamination_estimates_freebayes, \
        report_genotype_likelihood_max_freebayes, genotyping_max_iterations_freebayes, genotyping_max_banddepth_freebayes, posterior_integration_limits_freebayes, exclude_unobserved_genotypes_freebayes, \
        genotype_variant_threshold_freebayes, use_mapping_quality_freebayes, harmonic_indel_quality_freebayes, read_dependence_factor_freebayes, genotype_qualities_freebayes, debug_freebayes, dd_freebayes)

    emit:
        freebayes.out.collect()

}


workflow{
     variant_freebayes(channel.fromPath(params.bam), channel.fromPath(params.bai), 1)
}
