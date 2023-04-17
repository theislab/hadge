#!/usr/bin/env nextflow
'''
        path "sorted.bam"
        path "sorted.bam.bai"
        path "filtered.bam"
        path "filtered.bam.bai"
        path "no_dup.bam"
            sort_bam = samstool.out[0]
        sort_bamindex = samstool.out[1]
        
'''
nextflow.enable.dsl=2

process samstool{
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti/samtools", mode: 'copy'

    input:
        file bam

    output:
        path "samtools_${task.index}"
        

    script:
        """
        mkdir samtools_${task.index}
        samtools view -S -b -q 10 -F 3844 $bam > samtools_${task.index}/filtered.bam
        samtools index samtools_${task.index}/filtered.bam samtools_${task.index}/filtered.bam.bai
        umi_tools dedup --stdin=samtools_${task.index}/filtered.bam --extract-umi-method=tag --umi-tag=UR --cell-tag=CB --log=logfile > samtools_${task.index}/no_dup.bam
        samtools sort samtools_${task.index}/no_dup.bam -o samtools_${task.index}/sorted.bam
        samtools index samtools_${task.index}/sorted.bam samtools_${task.index}/sorted.bam.bai
        """
}


workflow data_preprocess{
    take:
        bam
    main:
        samstool(bam)
    emit:
        samstool.out

}

workflow{
    data_preprocess(channel.fromPath(params.bam))
}
