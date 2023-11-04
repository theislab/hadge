#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process vireo{
    publishDir "$projectDir/$params.outdir/$params.mode/gene_demulti/vireo", mode: 'copy'
    label 'big_mem'

    conda "aksarkar::vireosnp"

    input:
        each celldata
        each ndonor
        each donorfile
        each genoTag
        each noDoublet
        each nInit
        each extraDonor
        each extraDonorMode
        each forceLearnGT
        each ASEmode
        each noPlot
        each randSeed
        each cellRange
        each callAmbientRNAs
        each nproc
        each findVariant
        each vireo_out


    output:
        path "vireo_${task.index}"
    

    script:
        def cell_data = "-c $celldata"
        def n_donor =  ndonor != 'None'? "-N $ndonor" : ''
        def n_donor_yesno =  ndonor != 'None'? "$ndonor" : "Number of donors are not given"
        def donor = donorfile != 'None' ? "-d $donorfile" : ''
        def donor_data_name = donorfile != 'None' ? donorfile : 'Donor file is not given'
        def geno_tag = donorfile != 'None' ? "--genoTag $genoTag" : ''
        def no_doublet = noDoublet != 'False' ? "--noDoublet" : ''
        def n_init = "--nInit $nInit"
        def extra_donor = "--extraDonor $extraDonor"
        def extradonor_mode = extraDonorMode != 'distance' ? "--extraDonorMode $extraDonorMode" : ''
        def learnGT = (forceLearnGT != 'False' && donorfile != 'None')? "--forceLearnGT" : ''
        def learnGT_yesno = (forceLearnGT != 'False' && donorfile != 'None')? "$forceLearnGT" : 'False'
        def ase_mode = ASEmode != 'False' ? "--ASEmode" : ''
        def no_plot = noPlot != 'False' ? "--noPlot" : ''
        def random_seed = randSeed != 'None'? "--randSeed $randSeed" : ''
        def cell_range = cellRange != 'all'? "--cellRange $cellRange" : ''
        def call_ambient_rna = callAmbientRNAs != 'False' ? "--callAmbientRNAs" : ''
        def n_proc = "--nproc $nproc"

        """
        mkdir vireo_${task.index}
        mkdir vireo_${task.index}/${vireo_out}
        touch vireo_${task.index}/params.csv
        echo -e "Argument,Value \n cell_data,${celldata} \n n_donor,${n_donor_yesno} \n donor_data,${donor_data_name} \n genoTag,${genoTag} \n noDoublet,${noDoublet} \n nInit,${nInit} \n extraDonor,${extraDonor} \n extraDonorMode,${extraDonorMode} \n learnGT,${learnGT_yesno} \n ASEmode,${ASEmode} \n noPlot,${noPlot} \n randSeed,${randSeed} \n cellRange,${cellRange} \n callAmbientRNAs,${callAmbientRNAs} \n nproc,${nproc}" >> vireo_${task.index}/params.csv
        vireo ${cell_data} ${n_donor} $donor ${geno_tag} ${no_doublet} ${n_init} ${extra_donor} ${extradonor_mode} \
            $learnGT ${ase_mode} ${no_plot} ${random_seed} ${cell_range} ${call_ambient_rna} ${n_proc} \
            -o vireo_${task.index}/${vireo_out}
        if ([ "$donorfile" = "None" ]); then
            if ([ "$findVariant" = "True" ] || [ "$findVariant" = "vireo" ]); then
                GTbarcode -i vireo_${task.index}/${vireo_out}/GT_donors.vireo.vcf.gz -o vireo_${task.index}/${vireo_out}/filtered_variants.tsv ${randSeed}
            fi
        fi
        
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

workflow demultiplex_vireo{
    take:
        celldata
          
    main:
        ndonor = split_input(params.nsample)
        donorfile = split_input(params.vcf_donor)
        genoTag = split_input(params.genoTag)

        noDoublet = split_input(params.noDoublet)
        nInit = split_input(params.nInit)
        extraDonor = split_input(params.extraDonor)
        extraDonorMode = split_input(params.extraDonorMode)
        forceLearnGT = split_input(params.forceLearnGT)

        ASEmode = split_input(params.ASEmode)
        noPlot = split_input(params.noPlot)
        randSeed = split_input(params.randSeed)
        cellRange = split_input(params.cellRange)
        callAmbientRNAs = split_input(params.callAmbientRNAs)
        nproc = split_input(params.nproc)
        findVariant = split_input(params.findVariants)
        vireo_out = params.vireo_out
    
        vireo(celldata, ndonor, donorfile, genoTag, noDoublet, nInit, extraDonor, extraDonorMode, forceLearnGT, 
            ASEmode, noPlot, randSeed, cellRange, callAmbientRNAs, nproc, findVariant, vireo_out)
    

    emit:
        vireo.out.collect()
}
