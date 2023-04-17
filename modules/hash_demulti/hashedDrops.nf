#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process hashedDrops{
    publishDir "$projectDir/$params.outdir/$params.mode/hash_demulti/hashedDrops", mode:'copy'
    input:
        each raw_hto_matrix_dir
        each lower
        each niters
        each testAmbient
        each ignore
        each alpha
        each round
        each byRank
        each isCellFDR
        each objectOutEmptyDrops
        each assignmentOutEmptyDrops
	   
        each ambient
        each minProp
        each pseudoCount
        each constantAmbient
        each doubletNmads
        each doubletMin
        each doubletMixture
        each confidentNmads
        each confidenMin
        each combinations
        each objectOutHashedDrops
        each assignmentOutHashedDrops
    output:
        path "hashedDrops_${task.index}"
        
    script:
	    def testAmb = testAmbient != 'FALSE' ? " --testAmbient" : ''
        def rou = round != 'FALSE' ? " --round" : ''
        def constantAmb = constantAmbient != 'FALSE' ? " --constantAmbient" : ''
        def doubletMix = doubletMixture != 'FALSE' ?  " --doubletMixture" : ''

	    def ign = ignore != 'NULL' ? " --ignore ${ignore}" : ''
        def alp = alpha != 'NULL' ? " --alpha ${alpha}" : ''
        def byR = byRank != 'NULL' ? " --by.rank ${byRank}" : ''
        def amb = ambient != 'NULL' ? " --ambient ${ambient}" : ''
        def comb = combinations != 'NULL' ? " --combinations ${combinations}" : ''

        """
        mkdir hashedDrops_${task.index}
        dropletUtils.R --raw_hto_matrix_dir $raw_hto_matrix_dir --lower $lower --niters $niters --isCellFDR $isCellFDR --objectOutEmptyDrops $objectOutEmptyDrops --assignmentOutEmptyDrops $assignmentOutEmptyDrops --minProp $minProp --pseudoCount $pseudoCount --doubletNmads $doubletNmads --doubletMin $doubletMin --confidentNmads $confidentNmads --confidenMin $confidenMin --objectOutHashedDrops $objectOutHashedDrops --outputdir hashedDrops_${task.index} --assignmentOutHashedDrops ${assignmentOutHashedDrops}${testAmb}${ign}${alp}${rou}${byR}${constantAmb}${doubletMix}${amb}${comb}
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


workflow hashedDrops_hashing{
    main:
        raw_hto_matrix_dir = split_input(params.hto_matrix_hashedDrops)
        lower = split_input(params.lower)
        niters = split_input(params.niters)
        testAmbient = split_input(params.testAmbient)
        ignore = split_input(params.ignore_hashedDrops)
        alpha = split_input(params.alpha_hashedDrops)
        round = split_input(params.round)
        byRank = split_input(params.byRank)
        isCellFDR = split_input(params.isCellFDR)
        objectOutEmptyDrops = split_input(params.objectOutEmptyDrops)
        assignmentOutEmptyDrops = split_input(params.assignmentOutEmptyDrops)
       
        ambient = split_input(params.ambient)
        minProp = split_input(params.minProp)
        pseudoCount = split_input(params.pseudoCount)
        constantAmbient = split_input(params.constantAmbient)
        doubletNmads = split_input(params.doubletNmads)
        doubletMin = split_input(params.doubletMin)
        doubletMixture = split_input(params.doubletMixture)
        confidentNmads = split_input(params.confidentNmads)
        confidenMin = split_input(params.confidenMin)
        combinations = split_input(params.combinations)
        objectOutHashedDrops = split_input(params.objectOutHashedDrops)
        assignmentOutHashedDrops = split_input(params.assignmentOutHashedDrops)

        hashedDrops(raw_hto_matrix_dir, lower, niters, testAmbient, ignore, alpha, round, byRank, isCellFDR, objectOutEmptyDrops, assignmentOutEmptyDrops, ambient, minProp, pseudoCount, constantAmbient, doubletNmads, doubletMin, doubletMixture, confidentNmads, confidenMin, combinations, objectOutHashedDrops, assignmentOutHashedDrops)
        
    emit:
        hashedDrops.out.collect()
}

workflow{
    hashedDrops_hashing()
}
