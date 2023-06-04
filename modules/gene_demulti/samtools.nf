#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process samstool{
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode/gene_demulti/samtools", mode: 'copy'
    label 'big_mem'

    input:
        tuple val(sampleId), path(bam)

    output:
        path "samtools_${sampleId}"
        

    script:
        """
        mkdir samtools_${sampleId}
        samtools view -S -b -q 10 -F 3844 $bam > samtools_${sampleId}/filtered.bam
        samtools index samtools_${sampleId}/filtered.bam samtools_${sampleId}/filtered.bam.bai
        umi_tools dedup --stdin=samtools_${sampleId}/filtered.bam --extract-umi-method=tag --umi-tag=UR --cell-tag=CB --log=logfile > samtools_${sampleId}/no_dup.bam
        samtools sort samtools_${sampleId}/no_dup.bam -o samtools_${sampleId}/sorted.bam
        samtools index samtools_${sampleId}/sorted.bam samtools_${sampleId}/sorted.bam.bai
        """
}


workflow data_preprocess{
    take:
        input_list
    main:
        samstool(input_list)
    emit:
        samstool.out

}

