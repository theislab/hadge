#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process souporcell{
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode/gene_demulti/souporcell", mode: 'copy'
    label 'big_mem'

    container "shub://wheaton5/souporcell"
    
    input:
        tuple val(sampleId), path(bam), path(bam_index), path(barcodes), val(clusters), val(known_genotypes)
        path fasta
        val threads
        val ploidy
        val min_alt
        val min_ref
        val max_loci
        val restarts
        val common_variants
        val use_known_genotype
        val known_genotypes_sample_names
        val skip_remap
        val ignore
        val souporcell_out

    output:
        path "souporcell_${sampleId}"
    
    script:
        def bamfile = "-i $bam"
        def barcode = "-b $barcodes"
        def fastafile = "-f $fasta"
        def thread = "-t $threads"
        def cluster = "-k $clusters"
        def ploi = "--ploidy $ploidy"
        def minalt = "--min_alt ${min_alt}"
        def minref = "--min_ref ${min_ref}"
        def maxloci = "--max_loci ${max_loci}"
        def restart = restarts != 'None' ? "--restarts $restarts" : ''
        def commonvariant = (common_variants != 'None' & use_known_genotype != "True" & known_genotypes == 'None' )? "--common_variants ${common_variants}" : ''
        def commonvariant_name = (common_variants != 'None' & use_known_genotype != "True" & known_genotypes == 'None' ) ? common_variants : 'No common variants are given.'

        def knowngenotype = (known_genotypes != 'None' & use_known_genotype == "True") ? "--known_genotypes ${known_genotypes}" : ''
        def knowngenotype_name = (known_genotypes != 'None' & use_known_genotype == "True") ? known_genotypes : 'No known variants are given.'

        def knowngenotypes_sample = known_genotypes_sample_names != 'None' ? "--known_genotypes_sample_names ${known_genotypes_sample_names}" : ''
        def knowngenotype_sample_name = known_genotypes_sample_names != 'None' ? known_genotypes_sample_names : 'No known sample names are given.'

        def skipremap = skip_remap != 'False' ? "--skip_remap True" : ''
        def ign = ignore != 'False' ? "--ignore True" : ''
        def out = "souporcell_${sampleId}/${souporcell_out}"
        
        """
        mkdir souporcell_${sampleId}
        mkdir $out
        touch souporcell_${sampleId}/params.csv
        echo -e "Argument,Value \n bamfile,${bam} \n barcode,${barcodes} \n fasta,${fasta} \n threads,${threads} \n clusters,${clusters} \n ploidy,${ploidy} \n min_alt,${min_alt} \n min_ref,${min_ref} \n max_loci,${max_loci} \n restarts,${restarts} \n common_variant,${commonvariant_name} \n known_genotype,${knowngenotype_name} \n known_genotype_sample,${knowngenotype_sample_name} \n skip_remap,${skip_remap} \n ignore,${ignore} " >> souporcell_${sampleId}/params.csv
        
        souporcell_pipeline.py --threads ${task.cpus} $bamfile $barcode $fastafile $thread $cluster $ploi $minalt $minref $maxloci $restart \
            $commonvariant $knowngenotype $knowngenotypes_sample $skipremap $ign -o $out
        """
}

workflow demultiplex_souporcell{
    take:
        input_list
    main:
        fasta = params.fasta
        threads = params.threads
        ploidy = params.ploidy
        min_alt = params.min_alt
        min_ref = params.min_ref
        max_loci = params.max_loci
        restarts = params.restarts
        common_variants = params.common_variants_souporcell
        use_known_genotype = params.use_known_genotype
        known_genotypes_sample_names = params.known_genotypes_sample_names
        skip_remap = params.skip_remap
        ignore = params.ignore
        souporcell_out = params.souporcell_out
     
        souporcell(input_list, fasta, threads, ploidy, min_alt, min_ref, max_loci, restarts, common_variants, 
            use_known_genotype, known_genotypes_sample_names, skip_remap, ignore, souporcell_out)
            
    emit:
        souporcell.out.collect()
}
