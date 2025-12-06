/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_hadge_pipeline'

include { UNTAR as UNTAR_RNA                                       } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_HTO                                       } from '../modules/nf-core/untar/main'
include { RENAME_GENES_TO_FEATURES as RENAME_GENES_TO_FEATURES_RNA } from '../modules/local/rename_genes_to_features/main'
include { RENAME_GENES_TO_FEATURES as RENAME_GENES_TO_FEATURES_HTO } from '../modules/local/rename_genes_to_features/main'
include { EXTRACT_HASHES                                           } from '../modules/local/extract_hashes/main'

include { GENETIC_DEMULTIPLEXING     } from '../subworkflows/local/genetic_demultiplexing/main'
include { HASH_DEMULTIPLEXING        } from '../subworkflows/local/hash_demultiplexing/main'
include { CSVTK_JOIN as JOIN_RESULTS } from '../modules/nf-core/csvtk/join/main'
include { DONOR_MATCH                } from '../modules/local/donor_match/main'
include { FIND_VARIANTS              } from '../modules/local/find_variants/main'
include { SUBSET_GT_DONORS           } from '../modules/local/subset_gt_donors/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HADGE {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fasta // file: /path/to/genome.fasta

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    // ------------------------------ preprocessing start -------------------------------

    ch_rna = ch_samplesheet.map { meta, rna, _hto, _bam, _barcodes, _vcf -> [meta, rna] }
                    .branch { _meta, rna ->
                        tar: rna.endsWith('.tar.gz')
                        directory: true
                    }
    ch_hto = ch_samplesheet.map { meta, _rna, hto, _bam, _barcodes, _vcf -> [meta, hto] }
                    .branch { _meta, hto ->
                        tar: hto.endsWith('.tar.gz')
                        directory: true
                    }

    ch_remaining_input = ch_samplesheet.map { meta, _rna, _hto, bam, barcodes, vcf -> [meta, bam, barcodes, vcf] }

    UNTAR_RNA(ch_rna.tar)
    ch_versions = ch_versions.mix(UNTAR_RNA.out.versions)

    UNTAR_HTO(ch_hto.tar)
    ch_versions = ch_versions.mix(UNTAR_HTO.out.versions)

    ch_rna = ch_rna.directory.mix(UNTAR_RNA.out.untar)
    ch_hto = ch_hto.directory.mix(UNTAR_HTO.out.untar)

    // TODO remove completely if not used anymore
    // ch_rna = RENAME_GENES_TO_FEATURES_RNA(ch_rna)
    // ch_hto = RENAME_GENES_TO_FEATURES_HTO(ch_hto)

    // TODO maybe remove changes to extract hashes
    ch_hashes = EXTRACT_HASHES(ch_hto)

    ch_genetic = ch_samplesheet.map { meta, _rna, _hto, _bam, _barcodes, _vcf -> [meta] }
                        .join(ch_rna)
                        .join(ch_hto)
                        .join(ch_remaining_input)
                        .join(ch_hashes)
                        .map {meta, rna, hto, bam, barcodes, vcf, hashes -> [meta+[hashes: file(hashes).text.trim()], rna, hto, bam, barcodes, vcf] }

    ch_hashing = ch_genetic.map { meta, rna, hto, _bam, _barcodes, _vcf ->
        [meta, rna, hto]
    }

    ch_donor_match = ch_genetic.map { meta, _rna, _hto, _bam, barcodes, _vcf ->
        [meta, barcodes]
    }

    ch_find_variants = ch_donor_match.map { meta, barcodes -> [meta] }
    ch_subset_gt_donors = ch_donor_match.map { meta, barcodes -> [meta] }

    // ------------------------------- preprocessing end --------------------------------

    if (params.mode == 'genetic'){

        GENETIC_DEMULTIPLEXING(
            ch_genetic,
            params.genetic_tools.split(','),
            params.bam_qc,
            params.common_variants,
            fasta
        )

        ch_donor_match = ch_donor_match
            .join(GENETIC_DEMULTIPLEXING.out.summary_assignment)

        ch_versions = ch_versions.mix(GENETIC_DEMULTIPLEXING.out.versions)
    }
    else if (params.mode == 'hashing'){
        //TODO should mode hashing work with cell_genotype? nooo! -> maybe yes default would work
        // also maybe just for hashing

        HASH_DEMULTIPLEXING(
            ch_hashing,
            params.hash_tools.split(',')
        )

        ch_donor_match = ch_donor_match
            .join(HASH_DEMULTIPLEXING.out.summary_assignment)

        ch_versions = ch_versions.mix(HASH_DEMULTIPLEXING.out.versions)
    }
    else if ( params.mode == 'rescue' ){

        GENETIC_DEMULTIPLEXING(
            ch_genetic,
            params.genetic_tools.split(','),
            params.bam_qc,
            params.common_variants,
            fasta
        )

        HASH_DEMULTIPLEXING(
            ch_hashing,
            params.hash_tools.split(',')
        )

        JOIN_RESULTS(
            GENETIC_DEMULTIPLEXING.out.summary_assignment
                .join(HASH_DEMULTIPLEXING.out.summary_assignment)
                .map{meta, gene_summary, hash_summary ->
                    [meta, [gene_summary,hash_summary]]
                }
        )

        ch_donor_match = ch_donor_match
            .join(JOIN_RESULTS.out.csv)

        if ( params.find_variants ){
            ch_find_variants = ch_find_variants
                .join(GENETIC_DEMULTIPLEXING.out.gt_cells)
                .join(GENETIC_DEMULTIPLEXING.out.vireo_filtered_variants)
        }

        ch_versions = ch_versions.mix(GENETIC_DEMULTIPLEXING.out.versions)
        ch_versions = ch_versions.mix(HASH_DEMULTIPLEXING.out.versions)
        ch_versions = ch_versions.mix(JOIN_RESULTS.out.versions)
    }
    else if ( params.mode == 'donor_match' ){

        ['vireo_filtered_variants'].each { p ->
            if( !params[p] )
                error "Parameter '${p}' must be specified to run DONOR_MATCH with mode 'donor_match'"
            if( !file(params[p]).exists() )
                error "File specified for parameter '${p}' does not exist: ${params[p]}"
        }

        ch_donor_match = ch_donor_match.map{
            meta, barcodes ->
            [meta, barcodes, params.demultiplexing_result]
        }

        if ( params.find_variants ){

            ['cell_genotype', 'vireo_filtered_variants'].each { p ->
                if( !params[p] )
                    error "Parameter '${p}' must be specified to run FIND_VARIANTS with mode 'donor_match'"
                if( !file(params[p]).exists() )
                    error "File specified for parameter '${p}' does not exist: ${params[p]}"
            }

            ch_find_variants = ch_find_variants.map{ meta ->
                [meta, params.cell_genotype, params.vireo_filtered_variants]
            }
        }

    }

    if (params.match_donor) {

        DONOR_MATCH(ch_donor_match,
            params.match_donor_method1 ?: [],
            params.match_donor_method2 ?: []
        )

        // there only is a best_intersect_assignment_after_match output in donor_match and rescue mode
        if ( (params.mode == 'donor_match' | params.mode == 'rescue') && params.find_variants ){

            ch_find_variants = DONOR_MATCH.out.best_intersect_assignment_after_match
                .join(ch_find_variants)
                .join(ch_donor_match.map {
                        meta, barcode_whitelist, demultiplexing_result ->
                        [meta, demultiplexing_result]
                    }
                )

            FIND_VARIANTS(
                ch_find_variants,
                params.variant_count,
                params.variant_pct
            )

            // only vireo produces gt_donors
            if (params.genetic_tools && params.genetic_tools.split(',').contains('vireo')) {

                ch_subset_gt_donors = FIND_VARIANTS.out.donor_match_representative_variants
                    .map { meta, subset_variants ->
                        tuple(meta, subset_variants, 'donor_match')
                    }
                    .mix(
                        FIND_VARIANTS.out.vireo_representative_variants
                            .map { meta, subset_variants ->
                                tuple(meta, subset_variants, 'vireo')
                            }
                    )
                // ch_subset_gt_donors.view()

                ch_subset_gt_donors = ch_subset_gt_donors
                    .combine(GENETIC_DEMULTIPLEXING.out.gt_donors, by: 0)
                    .combine(DONOR_MATCH.out.best_donor_match, by: 0)
                    // .join(GENETIC_DEMULTIPLEXING.out.gt_donors)
                    // .join(DONOR_MATCH.out.best_donor_match)

                ch_subset_gt_donors.view()

                SUBSET_GT_DONORS(ch_subset_gt_donors)

            }

            ch_versions = ch_versions.mix(SUBSET_GT_DONORS.out.versions)
            ch_versions = ch_versions.mix(DONOR_MATCH.out.versions)
            ch_versions = ch_versions.mix(FIND_VARIANTS.out.versions)
        }


    }



    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'hadge_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
