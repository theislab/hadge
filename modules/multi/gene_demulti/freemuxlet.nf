#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process subset_bam_and_sort_vcf_based_on_reference{
    label 'small_mem'
    conda "-c conda-forge -c bioconda samtools=19.2 bedtools bcftools=1.19"
    tag "${sampleId}"
    input:
        tuple val(sampleId), path(sam), path(sam_index), path(group_list), val(nsample)
        path vcf

    output:
        tuple val(sampleId), path(sam), path(sam_index), path(group_list), val(nsample), path("${sampleId}__samples.sorted_as_in_bam.vcf"), emit: input
    
    script:
        """
            
            sort_vcf_same_as_bam.sh ${sam} ${vcf} > ${sampleId}__samples.sorted_as_in_bam.vcf        
        """
}

process freemuxlet {
    publishDir "$params.outdir/$sampleId/$params.mode/gene_demulti/freemuxlet", mode: 'copy'
    label 'small_mem'
    tag "${sampleId}"
    conda "bioconda::popscle"

    input:
        tuple val(sampleId), path(sam), path(sam_index), path(group_list), val(nsample), path(vcf)
        val tag_group
        val tag_UMI
        val sm
        val sm_list
        val sam_verbose
        val vcf_verbose
        val skip_umi
        val cap_BQ
        val min_BQ
        val min_MQ
        val min_TD
        val excl_flag
        val min_total
        val min_uniq
        val min_umi
        val min_snp
        val init_cluster
        val aux_files
        val verbose
        val doublet_prior
        val bf_thres
        val frac_init_clust
        val iter_init
        val keep_init_missing
        val freemuxlet_out

    output:
        path "freemuxlet_${sampleId}"

    script:
        def samfile = "--sam $sam"
        def taggroup = tag_group != 'None' ? "--tag-group ${tag_group}" : ''
        def tagUMI = tag_UMI != 'None' ? "--tag-UMI ${tag_UMI}" : ''
        def vcffile = "--vcf $vcf"
        def smlist = sm != 'None' ? "--sm $sm" : ''
        def sm_list_file = sm_list != 'None' ? "--sm-list ${sm_list}" : ''
        def sm_list_file_name = sm_list != 'None' ? sm_list : "No sm list file is given"
        def samverbose = "--sam-verbose ${sam_verbose}"
        def vcfverbose = "--vcf-verbose ${vcf_verbose}"
        def skipumi = skip_umi != "False" ? "--skip-umi" : ""
        def capBQ = "--cap-BQ ${cap_BQ}"
        def minBQ = "--min-BQ ${min_BQ}"
        def minMQ = "--min-MQ ${min_MQ}"
        def minTD = "--min-TD ${min_TD}"
        def exclflag = "--excl-flag ${excl_flag}"
        def grouplist =  "--group-list ${group_list}"
        def mintotal = "--min-total ${min_total}"
        def minuniq = "--min-uniq ${min_uniq}"
        def minumi = "--min-umi ${min_umi}"
        def minsnp = "--min-snp ${min_snp}"
        def initcluster = init_cluster != 'None' ? "--init-cluster ${init_cluster}" : ''
        def n_sample = "--nsample $nsample"
        def auxfiles = aux_files != 'False' ? "--aux-files" : ''
        def verbose_info = "--verbose $verbose"
        def doubletprior = "--doublet-prior ${doublet_prior}"
        def bfthres = "--bf-thres ${bf_thres}"
        def frac_init_cluster = "--frac-init-clust ${frac_init_clust}"
        def iterinit = "--iter-init ${iter_init}"
        def keepinit_missing = keep_init_missing != "False" ? "--keep-init-missing" : ''
        
        """
        mkdir freemuxlet_${sampleId}
        mkdir freemuxlet_${sampleId}/plp
        touch freemuxlet_${sampleId}/params.csv
        echo -e "Argument,Value \n samfile,${sam} \n tag_group,${tag_group} \n tag_UMI,${tag_UMI} \n vcf_file,${vcf} \n sm,${sm} \n sm_list_file,${sm_list_file_name} \n sam_verbose,${sam_verbose} \n vcf_verbose,${vcf_verbose} \n skip_umi,${skip_umi} \n cap_BQ,${cap_BQ} \n min_BQ,${min_BQ} \n min_MQ,${min_MQ} \n min_TD,${min_TD} \n excl_flag,${excl_flag} \n grouplist,${group_list} \n min_total,${min_total} \n min_uniq,${min_uniq} \n min_umi,${min_umi} \n min_snp,${min_snp} \n init_cluster,${init_cluster} \n nsample,${nsample} \n aux_files,${aux_files} \n verbose,${verbose} \n doublet_prior,${doublet_prior} \n bf_thres,${bf_thres} \n frac_init_clust,${frac_init_clust} \n iter_init,${iter_init} \n keep_init_missing,${keep_init_missing}" >> freemuxlet_${sampleId}/params.csv
        
        popscle dsc-pileup $samfile ${taggroup} ${tagUMI} $vcffile ${smlist} ${sm_list_file} ${samverbose} ${vcfverbose} \
            ${skipumi} ${capBQ} ${minBQ} ${minMQ} ${minTD} ${exclflag} ${grouplist} ${mintotal} ${minuniq} ${minsnp} \
            --out freemuxlet_${sampleId}/plp/${freemuxlet_out}
        popscle freemuxlet --plp freemuxlet_${sampleId}/plp/${freemuxlet_out} --out freemuxlet_${sampleId}/${freemuxlet_out} \
            ${initcluster} ${n_sample} ${auxfiles} ${verbose_info} ${doubletprior} ${bfthres} ${frac_init_cluster} ${iterinit} \
            ${keepinit_missing} ${capBQ} ${minBQ} ${grouplist} ${mintotal} ${minumi} ${minsnp}
        """
        
}



workflow demultiplex_freemuxlet{
    take:
        input_list
    main:
        vcf = params.common_variants_freemuxlet
        tag_group = params.tag_group
        tag_UMI = params.tag_UMI
        sm = params.sm
        sm_list = params.sm_list
        sam_verbose = params.sam_verbose
        vcf_verbose = params.vcf_verbose
        skip_umi = params.skip_umi
        cap_BQ = params.cap_BQ
        min_BQ = params.min_BQ
        min_MQ = params.min_MQ
        min_TD = params.min_TD
        excl_flag = params.excl_flag
        min_total = params.min_total
        min_uniq = params.min_uniq
        min_umi = params.min_umi
        min_snp = params.min_snp
        init_cluster = params.init_cluster
        aux_files = params.aux_files
        verbose = params.verbose
        doublet_prior = params.doublet_prior
        bf_thres = params.bf_thres
        frac_init_clust = params.frac_init_clust
        iter_init = params.iter_init
        keep_init_missing = params.keep_init_missing
        freemuxlet_out = params.freemuxlet_out



        subset_bam_and_sort_vcf_based_on_reference(input_list,vcf)
        input_list = subset_bam_and_sort_vcf_based_on_reference.out.input
        freemuxlet(input_list, tag_group, tag_UMI, sm, sm_list, sam_verbose, vcf_verbose, skip_umi, cap_BQ, min_BQ, min_MQ, 
            min_TD, excl_flag, min_total, min_uniq, min_umi, min_snp, init_cluster,aux_files, verbose, doublet_prior, bf_thres, 
            frac_init_clust, iter_init, keep_init_missing, freemuxlet_out)
    
    emit:
        freemuxlet.out.collect()
}
