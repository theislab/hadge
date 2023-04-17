#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process matchDonor{
    publishDir "$projectDir/$params.outdir/$params.mode/donor_match", mode: 'copy'
    input:
        path demultiplexing_result
        path barcode_whitelist
        val method1_name
        val method2_name
        val findVariants
        val cell_genotype
        val variant_count
        val variant_pct
        val vireo_parent_dir
    
    output:
        path "donor_match_${method1_name}_${method2_name}"

    script:
        def cell_genotype_path = ""
        if (findVariants == "True" | findVariants == "default"){
            if (cell_genotype == "False"){
                cell_genotype_path = "--cell_genotype $projectDir/$params.outdir/$params.mode/gene_demulti/cellSNP/cellsnp_1/*/cellSNP.cells.vcf.gz"
            }
            else{
                cell_genotype_path = "--cell_genotype $cell_genotype"
            }
        }
        def vireo_parent_path = ""
        if ( findVariants == 'vireo' | findVariants == 'True' ){
            if (params.mode == "donor_match" & vireo_parent_dir != 'False'){
                vireo_parent_path = "--vireo_parent_dir $vireo_parent_dir"
            }
            else{
                vireo_parent_path = "--vireo_parent_dir $projectDir/$params.outdir/$params.mode/gene_demulti/vireo/"
            }
            
        }
        def barcode_whitelist_path = ""
        if (barcode_whitelist != "False"){
          barcode_whitelist_path = "--barcode $barcode_whitelist"
        }
        """
        export R_MAX_VSIZE=100Gb
        outputdir=donor_match_${method1_name}_${method2_name}
        mkdir -p \$outputdir
        donor_match.R --result_csv $demultiplexing_result $barcode_whitelist_path --method1 $method1_name --method2 $method2_name  --findVariants $findVariants \
                $cell_genotype_path --variant_pct $variant_pct --variant_count $variant_count --outputdir \$outputdir $vireo_parent_path
                
        if ([ "$findVariants" != "False" ]); then
            best_method_vireo="\$(head -n 1 \$outputdir/best_method_vireo.txt)"
            if ([ "$params.mode" != "donor_match" ]); then
                donor_genotype="\$(find $projectDir/$params.outdir/$params.mode/gene_demulti/vireo/\$best_method_vireo -name GT_donors.vireo.vcf.gz | head -n 1)"
            else
                donor_genotype="\$(find $vireo_parent_dir/\$best_method_vireo -name GT_donors.vireo.vcf.gz | head -n 1)"
            fi
        
            if ([ "$findVariants" = "True" ] || [ "$findVariants" = "default" ]); then
                bcftools sort \$donor_genotype -Oz -o \$outputdir/compressed_sorted_donor_genotype.vcf.gz
                bcftools index \$outputdir/compressed_sorted_donor_genotype.vcf.gz
                bcftools filter \$outputdir/compressed_sorted_donor_genotype.vcf.gz -R \$outputdir/representative_variants_default.csv > \$outputdir/donor_genotype_subset_by_default.vcf
                bcftools reheader --samples \$outputdir/donor_match.csv -o \$outputdir/donor_genotype_subset_by_default_matched.vcf \$outputdir/donor_genotype_subset_by_default.vcf
                
            fi
            if ([ "$findVariants" = "True" ]); then
                bcftools filter \$outputdir/compressed_sorted_donor_genotype.vcf.gz -R \$outputdir/representative_variants_vireo.csv > \$outputdir/donor_genotype_subset_by_vireo.vcf
                bcftools reheader --samples \$outputdir/donor_match.csv -o \$outputdir/donor_genotype_subset_by_vireo_matched.vcf \$outputdir/donor_genotype_subset_by_vireo.vcf
            fi
            if ([ "$findVariants" = "vireo" ]); then
                bcftools sort \$donor_genotype -Oz -o \$outputdir/compressed_sorted_donor_genotype.vcf.gz
                bcftools index \$outputdir/compressed_sorted_donor_genotype.vcf.gz
                bcftools filter \$outputdir/compressed_sorted_donor_genotype.vcf.gz -R \$outputdir/representative_variants_vireo.csv > \$outputdir/donor_genotype_subset_by_vireo.vcf
                bcftools reheader --samples \$outputdir/donor_match.csv -o \$outputdir/donor_genotype_subset_by_vireo_matched.vcf \$outputdir/donor_genotype_subset_by_vireo.vcf
            fi
        fi


        """
}

                
workflow donor_match{
    take:
        demultiplexing_result
    main:
        matchDonor(demultiplexing_result, params.barcodes, params.match_donor_method1, params.match_donor_method2, params.findVariants, params.celldata,
            params.variant_count, params.variant_pct, params.vireo_parent_dir)
}
