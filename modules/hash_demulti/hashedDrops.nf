#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process hashedDrops{
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode/hash_demulti/hashedDrops", mode:'copy'
    label 'small_mem'

    conda "conda-forge::r-seurat conda-forge::r-argparse bioconda::bioconductor-dropletutils"
    
    input:
        tuple val(sampleId), path(raw_hto_matrix_dir)
        val lower
        val niters
        val testAmbient
        val ignore
        val alpha
        val round
        val byRank
        val isCellFDR
        val objectOutEmptyDrops
        val assignmentOutEmptyDrops
	   
        val ambient
        val minProp
        val pseudoCount
        val constantAmbient
        val doubletNmads
        val doubletMin
        val doubletMixture
        val confidentNmads
        val confidenMin
        val combinations
        val objectOutHashedDrops
        val assignmentOutHashedDrops
    output:
        path "hashedDrops_${sampleId}"
        
    script:
	    def testAmb = testAmbient != 'False' ? " --testAmbient" : ''
        def rou = round != 'False' ? " --round" : ''
        def constantAmb = constantAmbient != 'False' ? " --constantAmbient" : ''
        def doubletMix = doubletMixture != 'False' ?  " --doubletMixture" : ''
	    def ign = ignore != "None" ? " --ignore ${ignore}" : ''
        def alp = alpha != "None" ? " --alpha ${alpha}" : ''
        def byR = byRank != "None" ? " --by.rank ${byRank}" : ''
        def amb = ambient != 'False' ? " --ambient" : ''
        def comb = combinations != "None" ? " --combinations ${combinations}" : ''

        """
        mkdir hashedDrops_${sampleId}
        dropletUtils.R --raw_hto_matrix_dir $raw_hto_matrix_dir --lower $lower --niters $niters \
            --isCellFDR $isCellFDR --objectOutEmptyDrops $objectOutEmptyDrops \
            --assignmentOutEmptyDrops $assignmentOutEmptyDrops --minProp $minProp --pseudoCount $pseudoCount \
            --doubletNmads $doubletNmads --doubletMin $doubletMin --confidentNmads $confidentNmads \
            --confidenMin $confidenMin --objectOutHashedDrops $objectOutHashedDrops \
            --outputdir hashedDrops_${sampleId} --assignmentOutHashedDrops ${assignmentOutHashedDrops} \
            ${testAmb}${ign}${alp}${rou}${byR}${constantAmb}${doubletMix}${amb}${comb}
        """

}

workflow hashedDrops_hashing{
    take:
        input_list
    main:
        lower = params.lower
        niters = params.niters
        testAmbient = params.testAmbient
        ignore = params.ignore_hashedDrops
        alpha = params.alpha_hashedDrops
        round = params.round
        byRank = params.byRank
        isCellFDR = params.isCellFDR
        objectOutEmptyDrops = params.objectOutEmptyDrops
        assignmentOutEmptyDrops = params.assignmentOutEmptyDrops
       
        ambient = params.ambient
        minProp = params.minProp
        pseudoCount = params.pseudoCount
        constantAmbient = params.constantAmbient
        doubletNmads = params.doubletNmads
        doubletMin = params.doubletMin
        doubletMixture = params.doubletMixture
        confidentNmads = params.confidentNmads
        confidenMin = params.confidenMin
        combinations = params.combinations
        objectOutHashedDrops = params.objectOutHashedDrops
        assignmentOutHashedDrops = params.assignmentOutHashedDrops

        hashedDrops(input_list, lower, niters, testAmbient, ignore, alpha, round, byRank, isCellFDR, 
            objectOutEmptyDrops, assignmentOutEmptyDrops, ambient, minProp, pseudoCount, constantAmbient, 
            doubletNmads, doubletMin, doubletMixture, confidentNmads, confidenMin, combinations, objectOutHashedDrops, 
            assignmentOutHashedDrops)
        
    emit:
        hashedDrops.out.collect()
}

workflow{
    hashedDrops_hashing()
}
