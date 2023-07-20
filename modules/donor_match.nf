#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process matchDonor{
    publishDir "$projectDir/$params.outdir/$sampleId/$params.mode", mode: 'copy'
    label 'big_mem'
    
    input:
        tuple val (sampleId), val(ndonor), path(barcode_whitelist), val(cell_genotype), val(vireo_parent_dir), path(demultiplexing_result)
        val method1_name
        val method2_name
        val findVariants
        val variant_count
        val variant_pct
    
    output:
        tuple val(sampleId), path("donor_match_${sampleId}")

    script:
        def two_method = (method1_name != "None" & method2_name != "None") ? "--method1 $method1_name --method2 $method2_name" : ""
        
        def cell_genotype_path = ""
        if (findVariants == "True" | findVariants == "default"){
            cell_genotype_path = cell_genotype != "None" ? "--cell_genotype $cell_genotype" : \
                "--cell_genotype $projectDir/$params.outdir/$sampleId/$params.mode/gene_demulti/cellSNP/cellsnp_1/*/cellSNP.cells.vcf.gz"
        }

        def vireo_parent_path = ""
        if ( findVariants == 'vireo' | findVariants == 'True' ){
            vireo_parent_path = (params.mode == "donor_match" & vireo_parent_dir != 'None') ? "--vireo_parent_dir $vireo_parent_dir" : "--vireo_parent_dir $projectDir/$params.outdir/$sampleId/$params.mode/gene_demulti/vireo/"
        }
        def barcode_whitelist_path = "--barcode $barcode_whitelist"
        """
        export R_MAX_VSIZE=100Gb
        outputdir=donor_match_${sampleId}
        mkdir -p \$outputdir
        donor_match.R --result_csv $demultiplexing_result $barcode_whitelist_path --findVariants $findVariants \
                $cell_genotype_path --variant_pct $variant_pct --variant_count $variant_count --ndonor $ndonor \
                $two_method --outputdir \$outputdir $vireo_parent_path
                
        if ([ "$findVariants" != "False" ]); then
            best_method_vireo="\$(head -n 1 \$outputdir/best_method_vireo.txt)"
            if ([ "$params.mode" != "donor_match" ]); then
                donor_genotype="\$(find $projectDir/$params.outdir/$sampleId/$params.mode/gene_demulti/vireo/\$best_method_vireo -name GT_donors.vireo.vcf.gz | head -n 1)"
            else
                donor_genotype="\$(find $vireo_parent_dir/\$best_method_vireo -name GT_donors.vireo.vcf.gz | head -n 1)"
            fi
        
            if ([ "$findVariants" = "True" ] || [ "$findVariants" = "default" ]); then
                gunzip -c \$donor_genotype > \$outputdir/GT_donors.vireo.vcf
                if ([ \$(grep "^##config" \$outputdir/GT_donors.vireo.vcf | wc -l) == 0 ]); then
                    bcftools view --header-only $cell_genotype | grep "^##" | grep -v "^##bcftools" > \$outputdir/donor_with_header.vcf
                    grep -v "^##" \$outputdir/GT_donors.vireo.vcf >> \$outputdir/donor_with_header.vcf
                    bcftools sort \$outputdir/donor_with_header.vcf -Oz -o \$outputdir/compressed_sorted_donor_genotype.vcf.gz
                    rm \$outputdir/donor_with_header.vcf
                else
                    bcftools sort \$donor_genotype -Oz -o \$outputdir/compressed_sorted_donor_genotype.vcf.gz
                fi
                bcftools index \$outputdir/compressed_sorted_donor_genotype.vcf.gz
                bcftools filter \$outputdir/compressed_sorted_donor_genotype.vcf.gz -R \$outputdir/donor_specific_variants.csv > \$outputdir/donor_genotype_subset_by_default.vcf
                bcftools reheader --samples \$outputdir/donor_match.csv -o \$outputdir/donor_genotype_subset_by_default_matched.vcf \$outputdir/donor_genotype_subset_by_default.vcf

                rm \$outputdir/GT_donors.vireo.vcf
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
            rm \$outputdir/best_method_vireo.txt
            rm \$outputdir/compressed_sorted_donor_genotype.vcf.gz
            rm \$outputdir/compressed_sorted_donor_genotype.vcf.gz.csi
        fi
        
        """
}

                
workflow donor_match{
    take:
        demultiplexing_result
    main:
        matchDonor(demultiplexing_result, params.match_donor_method1, params.match_donor_method2, 
            params.findVariants, params.variant_count, params.variant_pct)
    emit:
        matchDonor.out
}
