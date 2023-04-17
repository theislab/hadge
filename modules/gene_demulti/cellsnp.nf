#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process cellSNP{
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti/cellSNP", mode: 'copy'
    
    input:
        val samFile_cellSNP
        file samFileList_cellSNP
        file regionsVCF_cellSNP
        file targetsVCF_cellSNP
        each barcodeFile_cellSNP
        file sampleList_cellSNP
        val sampleIDs_cellSNP
        val genotype_cellSNP
        val gzip_cellSNP
        val printSkipSNPs_cellSNP
        val nproc_cellSNP
        file refseq_cellSNP
        val chrom_cellSNP
        val cellTAG_cellSNP
        val UMItag_cellSNP
        val minCOUNT_cellSNP
        val minMAF_cellSNP
        val doubletGL_cellSNP
        val inclFLAG_cellSNP
        val exclFLAG_cellSNP
        val minLEN_cellSNP
        val minMAPQ_cellSNP
        val countORPHAN_cellSNP
        val cellsnp_out


    output:
        path "cellsnp_${task.index}"

    script:
        def samFile = samFile_cellSNP != 'no_samFile' ? "--samFile ${samFile_cellSNP}" : ''
        def samFileList = samFileList_cellSNP.name != 'no_samFileList' ? "--samFileList ${samFileList_cellSNP}" : ''
        def regionsVCF = regionsVCF_cellSNP.name != 'no_regionsVCF' ? "--regionsVCF ${regionsVCF_cellSNP}" : ''
        def targetsVCF = targetsVCF_cellSNP.name != 'no_targetsVCF' ? "--targetsVCF ${targetsVCF_cellSNP}" : ''
        def barcodeFile = barcodeFile_cellSNP != 'no_barcodeFile' ? "--barcodeFile ${barcodeFile_cellSNP}" : ''
        def sampleList = sampleList_cellSNP.name != 'no_sampleList' ? "--sampleList ${sampleList_cellSNP}" : ''
        def sampleIDs = sampleIDs_cellSNP != 'no_sampleIDs' ? "--sampleIDs ${sampleIDs_cellSNP}" : ''
        def genotype = genotype_cellSNP != 'False' ? "--genotype" : ''
        def gzip = gzip_cellSNP != 'False' ? "--gzip" : ''
        def printSkipSNPs = printSkipSNPs_cellSNP != 'False' ? "--printSkipSNPs" : ''
        def nproc = nproc_cellSNP != 'False' ? "--nproc ${nproc_cellSNP}" : ''
        def refseq = refseq_cellSNP.name != 'no_refseq' ? "--refseq ${refseq_cellSNP}" : ''
        def chrom = chrom_cellSNP != 'False' ? "--chrom ${chrom_cellSNP}" : ''
        def cellTAG = "--cellTAG ${cellTAG_cellSNP}"
        def UMItag = "--UMItag ${UMItag_cellSNP}"
        def minCOUNT = "--minCOUNT ${minCOUNT_cellSNP}"
        def minMAF = "--minMAF ${minMAF_cellSNP}"
        def doubletGL = doubletGL_cellSNP != 'False' ? "--doubletGL" : ''
        def inclFLAG = inclFLAG_cellSNP != 'False' ? "--inclFLAG ${inclFLAG_cellSNP}" : ''
        def exclFLAG = exclFLAG_cellSNP != 'False' ? "--exclFLAG ${exclFLAG_cellSNP}" : ''
        def minLEN = "--minLEN ${minLEN_cellSNP}"
        def minMAPQ = "--minMAPQ ${minMAPQ_cellSNP}"
        def countORPHAN = countORPHAN_cellSNP != 'False' ? "--countORPHAN" : ''
        def out = "cellsnp_${task.index}/${cellsnp_out}"

        """
        mkdir cellsnp_${task.index}
        mkdir $out
        touch cellsnp_${task.index}/params.csv
        echo -e "Argument,Value \n samfile,${samFile_cellSNP} \n samFileList,${samFileList_cellSNP} \n regionsVCF, ${regionsVCF_cellSNP} \n targetsVCF, ${targetsVCF_cellSNP} \n barcodeFile, ${barcodeFile_cellSNP} \n sampleList, ${sampleList_cellSNP} \n sampleIDs, ${sampleIDs_cellSNP} \n genotype, ${genotype_cellSNP} \n gzip, ${gzip_cellSNP} \n printSkipSNPs, ${printSkipSNPs_cellSNP} \n nproc, ${nproc_cellSNP} \n refseq, ${refseq_cellSNP} \n chrom, ${chrom_cellSNP} \n cellTAG, ${cellTAG_cellSNP} \n UMItag, ${UMItag_cellSNP} \n minCOUNT, ${minCOUNT_cellSNP} \n minMAF, ${minMAF_cellSNP} \n doubletGL, ${doubletGL_cellSNP} \n inclFLAG, ${inclFLAG_cellSNP} \n exclFLAG, ${exclFLAG_cellSNP} \n minLEN, ${minLEN_cellSNP} \n minMAPQ, ${minMAPQ_cellSNP} \n countORPHAN, ${countORPHAN_cellSNP}" >> cellsnp_${task.index}/params.csv
        cellsnp-lite $samFile $samFileList $regionsVCF $targetsVCF $barcodeFile $sampleList $sampleIDs $genotype $gzip $printSkipSNPs $nproc $refseq $chrom $cellTAG $UMItag $minCOUNT $minMAF $doubletGL $inclFLAG $exclFLAG $minLEN $minMAPQ $countORPHAN --outDir $out
        cd $out
        gunzip -c cellSNP.cells.vcf.gz > cellSNP.cells.vcf
        """
}


def split_input(input){
    if (input =~ /;/ ){
        Channel.from(input).map{ return it.tokenize(';')}.flatten()
    }
    else{
        Channel.from(input)
    }
}

workflow variant_cellSNP{
    take:
        samFile
        barcodeFile
    main:
        samFileList = channel.fromPath(params.samFileList)
        regionsVCF = channel.fromPath(params.regionsVCF)
        targetsVCF =  channel.fromPath(params.targetsVCF)
        sampleList = channel.fromPath(params.sampleList)
        sampleIDs = channel.value(params.sampleIDs)
        genotype_cellSNP = channel.value(params.genotype_cellSNP)
        gzip_cellSNP = channel.value(params.gzip_cellSNP)
        printSkipSNPs = channel.value(params.printSkipSNPs)
        nproc_cellSNP = channel.value(params.nproc_cellSNP)
        refseq_cellSNP = channel.fromPath(params.refseq_cellSNP)
        chrom = channel.value(params.chrom)
        cellTAG = channel.value(params.cellTAG)
        UMItag = channel.value(params.UMItag)
        minCOUNT = channel.value(params.minCOUNT)
        minMAF = channel.value(params.minMAF)
        doubletGL = channel.value(params.doubletGL)
        inclFLAG = channel.value(params.inclFLAG)
        exclFLAG = channel.value(params.exclFLAG)
        minLEN = channel.value(params.minLEN)
        minMAPQ = channel.value(params.minMAPQ)
        countORPHAN = channel.value(params.countORPHAN)
        cellsnp_out = channel.value(params.cellsnp_out)
        cellSNP(samFile, samFileList, regionsVCF, targetsVCF, barcodeFile, sampleList, sampleIDs, genotype_cellSNP, gzip_cellSNP, printSkipSNPs, nproc_cellSNP, refseq_cellSNP, chrom, cellTAG, UMItag, minCOUNT, minMAF, doubletGL, inclFLAG, exclFLAG, minLEN, minMAPQ, countORPHAN, cellsnp_out)
    emit:
        cellSNP.out
}

workflow{
     variant_cellSNP(channel.fromPath(params.bam))
}
