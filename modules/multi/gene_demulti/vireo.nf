#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process vireo{
    publishDir "$params.outdir/$sampleId/$params.mode/gene_demulti/vireo", mode: 'copy'
    label 'big_mem'
    tag "${sampleId}"
    conda "aksarkar::vireosnp"

    input:
        tuple val(sampleId), path(celldata), val(ndonor), val(donorfile)
        val genoTag
        val noDoublet
        val nInit
        val extraDonor
        val extraDonorMode
        val forceLearnGT
        val ASEmode
        val noPlot
        val randSeed
        val cellRange
        val callAmbientRNAs
        val nproc
        val findVariant
        val vireo_out


    output:
        path "vireo_${sampleId}"
    

    script:
        def cell_data = "-c $celldata"
        def celldata_name = celldata.baseName
        def n_donor =  ndonor != 'None'? "-N $ndonor" : ''
        def n_donor_yesno =  ndonor != 'None'? "$ndonor" : "Number of donors are not given"
        def donor = donorfile != 'None' ? "-d no_prefix_chr.vcf" : ''
        def donor_preprocess = donorfile != 'None' ? "bcftools view $donorfile | awk '{gsub(/^chr/,\"\"); print}' | awk '{gsub(/ID=chr/,\"ID=\"); print}' > no_prefix_chr.vcf" : ''
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
        ${donor_preprocess}
        mkdir vireo_${sampleId}
        mkdir vireo_${sampleId}/${vireo_out}
        touch vireo_${sampleId}/params.csv
        echo -e "Argument,Value \n cell_data,${celldata_name} \n n_donor,${n_donor_yesno} \n donor_data,${donor_data_name} \n genoTag,${genoTag} \n noDoublet,${noDoublet} \n nInit,${nInit} \n extraDonor,${extraDonor} \n extraDonorMode,${extraDonorMode} \n learnGT,${learnGT_yesno} \n ASEmode,${ASEmode} \n noPlot,${noPlot} \n randSeed,${randSeed} \n cellRange,${cellRange} \n callAmbientRNAs,${callAmbientRNAs} \n nproc,${nproc}" >> vireo_${sampleId}/params.csv
        
        vireo ${cell_data} ${n_donor} $donor ${geno_tag} ${no_doublet} ${n_init} ${extra_donor} ${extradonor_mode} $learnGT \
            ${ase_mode} ${no_plot} ${random_seed} ${cell_range} ${call_ambient_rna} ${n_proc} -o vireo_${sampleId}/${vireo_out}
        
        if ([ "$donor_data_name" = "Donor file is not given" ]); then
            if ([ "$findVariant" = "True" ] || [ "$findVariant" = "vireo" ]); then
                GTbarcode -i vireo_${sampleId}/${vireo_out}/GT_donors.vireo.vcf.gz -o vireo_${sampleId}/${vireo_out}/filtered_variants.tsv ${randSeed}
            fi
        fi
        
        """

}

workflow demultiplex_vireo{
    take:
        input_list
          
    main:
        genoTag = params.genoTag
        noDoublet = params.noDoublet
        nInit = params.nInit
        extraDonor = params.extraDonor
        extraDonorMode = params.extraDonorMode
        forceLearnGT = params.forceLearnGT
        ASEmode = params.ASEmode
        noPlot = params.noPlot
        randSeed = params.randSeed
        cellRange = params.cellRange
        callAmbientRNAs = params.callAmbientRNAs
        nproc = params.nproc
        findVariant = params.findVariants
        vireo_out = params.vireo_out
    
        vireo(input_list, genoTag, noDoublet, nInit, extraDonor, extraDonorMode, forceLearnGT, ASEmode, noPlot, randSeed, \
            cellRange, callAmbientRNAs, nproc, findVariant, vireo_out)
    

    emit:
        vireo.out.collect()
}
