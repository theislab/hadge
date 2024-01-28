#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
process bcftools {
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti/bcftools", mode: 'copy'
    label 'big_mem'

    conda 'bioconda::bcftools=1.9'

    input:
        val vcf
    output:
        path "bcftools_${task.index}"

    script:
        vcf_files = vcf.join(' ')
        """
        mkdir bcftools_${task.index}
        bcftools concat -o bcftools_${task.index}/total_chroms.vcf ${vcf_files}
        cd bcftools_${task.index}
        bcftools sort total_chroms.vcf -o sorted_total_chroms.vcf
        bcftools filter -i '%QUAL > 30' sorted_total_chroms.vcf -o filtered_sorted_total_chroms.vcf
        """
}
workflow filter_variant {
    take:
        vcf
    main:
        bcftools(vcf)
    emit:
        bcftools.out
}

