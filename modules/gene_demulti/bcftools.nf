#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process bcftools{
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti/bcftools", mode: 'copy'
    label 'big_mem'
    
    input:
        val vcf
        val concat
        val filter
    output:
        path "bcftools_${task.index}"

    script:
        vcf_files = ""
        if (concat == "True"){
            for(v : vcf) {
                vcf_files = v + " " + vcf_files
            }
            """
            mkdir bcftools_${task.index}
            bcftools concat -o bcftools_${task.index}/total_chroms.vcf ${vcf_files}
            cd bcftools_${task.index}
            bcftools sort total_chroms.vcf -o sorted_total_chroms.vcf
            if [[ "$filter" == "True" ]]
            then
                bcftools filter -i '%QUAL > 30' sorted_total_chroms.vcf -o filtered_sorted_total_chroms.vcf
            fi
            """
        }
        else if (filter == "True"){
            """
            mkdir bcftools_${task.index}
            bcftools filter -i '%QUAL > 30' $vcf -o bcftools_${task.index}/filtered.vcf
            """
        }
}
workflow filter_variant{
    take:
        vcf
        concat
        filter
    main:
        bcftools(vcf, concat, filter)
    emit:
        bcftools.out
}

