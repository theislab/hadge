#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process bcftools{
    publishDir "$params.outdir/$sampleId/$params.mode/gene_demulti/bcftools", mode: 'copy'
    label 'big_mem'
    tag "${sampleId}"
    conda "bioconda::bcftools=1.9"
    
    input:
        tuple val(sampleId), val(vcf_list)
    output:
        tuple val(sampleId), path("bcftools_${sampleId}")

    script:
        vcf_files = vcf_list.join(" ")
        """
        mkdir bcftools_${sampleId}
        bcftools concat -o bcftools_${sampleId}/total_chroms.vcf ${vcf_files}
        cd bcftools_${sampleId}
        bcftools sort total_chroms.vcf -o sorted_total_chroms.vcf
        bcftools filter -i '%QUAL > 30' sorted_total_chroms.vcf -o filtered_sorted_total_chroms.vcf
        """
       
}
workflow filter_variant{
    take:
        input_list
    main:
        bcftools(input_list)
    emit:
        bcftools.out
}

