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
include { EXTRACT_HASHES                                           } from '../modules/local/extract_hashes/main'

include { GENETIC_DEMULTIPLEXING                    } from '../subworkflows/local/genetic_demultiplexing/main'
include { HASH_DEMULTIPLEXING                       } from '../subworkflows/local/hash_demultiplexing/main'
include { CREATE_ANNDATA_MUDATA                     } from '../modules/local/create_anndata_mudata/main'
include { CSVTK_JOIN as JOIN_RESULTS_ASSIGNMENT     } from '../modules/nf-core/csvtk/join/main'
include { CSVTK_JOIN as JOIN_RESULTS_CLASSIFICATION } from '../modules/nf-core/csvtk/join/main'
include { DONOR_MATCH                               } from '../modules/local/donor_match/main'
include { FIND_VARIANTS                             } from '../modules/local/find_variants/main'
include { SUBSET_GT_DONORS                          } from '../modules/local/subset_gt_donors/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def checkParams(String paramName, String process, String mode, boolean isFile) {
    def value = params[paramName]

    if( !value )
        error "Parameter '${paramName}' must be specified to run ${process} with mode '${mode}'"

    if( !value && !mode )
        error "Parameter '${paramName}' must be specified to run ${process}"

    if( isFile && !file(value).exists() )
        error "File specified for parameter '${paramName}' does not exist: ${value}"

    return true
}

workflow HADGE {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fasta // file: /path/to/genome.fasta

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // ------------------------------ preprocessing start -------------------------------
    if ( params.mode != 'donor_match' ){
        ch_rna = ch_samplesheet.map { meta, rna, _hto, _bam, _barcodes, _vcf -> [meta, rna] }
                        .branch { _meta, rna ->
                            tar: rna.endsWith('.tar.gz')
                            directory: true
                        }


        ch_hto = ch_samplesheet.map { meta, _rna, hto, _bam, _barcodes, _vcf -> [meta, hto] }
                        .branch { _meta, hto ->
                            tar: hto != null && hto.endsWith('.tar.gz')
                            directory: true
                        }

        ch_remaining_input = ch_samplesheet.map { meta, _rna, _hto, bam, barcodes, vcf -> [meta, bam, barcodes, vcf] }

        // @nictru do I have to track versions of both modules even tough it is from the same module?
        UNTAR_RNA(ch_rna.tar)
        ch_versions = ch_versions.mix(UNTAR_RNA.out.versions)

        UNTAR_HTO(ch_hto.tar)
        ch_versions = ch_versions.mix(UNTAR_HTO.out.versions)

        ch_rna = ch_rna.directory.mix(UNTAR_RNA.out.untar)
        ch_hto = ch_hto.directory.mix(UNTAR_HTO.out.untar)

        // hto can be null in genetic mode
        ch_hashes_non_null = EXTRACT_HASHES(ch_hto.filter { _meta, hto -> hto != null })
        ch_hashes_null = ch_hto.filter { _meta, hto -> hto == null }
        ch_hashes = ch_hashes_non_null.mix(ch_hashes_null)

        ch_preprocessed = ch_samplesheet.map { meta, _rna, _hto, _bam, _barcodes, _vcf -> [meta] }
                            .join(ch_rna)
                            .join(ch_hto)
                            .join(ch_remaining_input)
                            .join(ch_hashes)
                            .map {meta, rna, hto, bam, barcodes, vcf, hashes ->
                            if(hashes != null && meta.hto_names == []){meta += [hto_names: file(hashes).text.trim()]}
                            [meta, rna, hto, bam, barcodes, vcf]
                            }

        // create channels for deconvolution tools
        ch_genetic = ch_preprocessed.map { meta, rna, _hto, bam, barcodes, vcf ->
            [meta, rna, bam, barcodes, vcf]
        }

        ch_hashing = ch_preprocessed.map { meta, rna, hto, _bam, _barcodes, _vcf ->
            [meta, rna, hto]
        }

        ch_create_anndata_mudata = ch_preprocessed.map { meta, rna, hto, _bam, _barcodes, _vcf -> [meta, rna, hto] }
    }else{
        // meta changes when extracting hashes
        ch_preprocessed = ch_samplesheet
    }

    // channels for donor matching
    ch_donor_match = ch_preprocessed.map { meta, _rna, _hto, _bam, barcodes, _vcf ->
        [meta, barcodes]
    }
    ch_find_variants = ch_donor_match.map { meta, _barcodes -> [meta] }
    ch_subset_gt_donors = ch_donor_match.map { meta, _barcodes -> [meta] }

    // ------------------------------- preprocessing end --------------------------------

    if (params.mode == 'genetic'){

        GENETIC_DEMULTIPLEXING(
            ch_genetic,
            params.genetic_tools.split(','),
            params.bam_qc,
            params.common_variants,
            fasta
        )

        ch_create_anndata_mudata = ch_create_anndata_mudata
            .join(GENETIC_DEMULTIPLEXING.out.summary_assignment)
            .join(GENETIC_DEMULTIPLEXING.out.summary_classification)
            .map { meta, rna, hto, gene_a, gene_c ->
                [meta, rna, hto, gene_a, gene_c, [], []]
            }

        ch_donor_match = ch_donor_match
            .join(GENETIC_DEMULTIPLEXING.out.summary_assignment)

        ch_versions = ch_versions.mix(GENETIC_DEMULTIPLEXING.out.versions)
    }

    else if (params.mode == 'hashing'){

        HASH_DEMULTIPLEXING(
            ch_hashing,
            params.hash_tools.split(',')
        )

        ch_create_anndata_mudata = ch_create_anndata_mudata
            .map { meta, rna, hto -> [meta, rna, hto, [], []]}
            .join(HASH_DEMULTIPLEXING.out.summary_assignment)
            .join(HASH_DEMULTIPLEXING.out.summary_classification)

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

        JOIN_RESULTS_ASSIGNMENT(
            GENETIC_DEMULTIPLEXING.out.summary_assignment
                .join(HASH_DEMULTIPLEXING.out.summary_assignment)
                .map{meta, gene_summary, hash_summary ->
                    [meta, [gene_summary,hash_summary]]
                }
        )

        JOIN_RESULTS_CLASSIFICATION(
            GENETIC_DEMULTIPLEXING.out.summary_classification
                .join(HASH_DEMULTIPLEXING.out.summary_classification)
                .map{meta, gene_summary, hash_summary ->
                    [meta, [gene_summary,hash_summary]]
                }
        )

        ch_create_anndata_mudata = ch_create_anndata_mudata
            .join(GENETIC_DEMULTIPLEXING.out.summary_assignment)
            .join(GENETIC_DEMULTIPLEXING.out.summary_classification)
            .join(HASH_DEMULTIPLEXING.out.summary_assignment)
            .join(HASH_DEMULTIPLEXING.out.summary_classification)

        ch_donor_match = ch_donor_match
            .join(JOIN_RESULTS_ASSIGNMENT.out.csv)

        if ( params.find_variants ){
            ch_find_variants = ch_find_variants
                .join(GENETIC_DEMULTIPLEXING.out.gt_cells)
                .join(GENETIC_DEMULTIPLEXING.out.vireo_filtered_variants)
        }

        ch_versions = ch_versions.mix(GENETIC_DEMULTIPLEXING.out.versions)
        ch_versions = ch_versions.mix(HASH_DEMULTIPLEXING.out.versions)
        // @nictru do I have to track versions of both modules even tough it is from the same module?
        ch_versions = ch_versions.mix(JOIN_RESULTS_ASSIGNMENT.out.versions)
        ch_versions = ch_versions.mix(JOIN_RESULTS_CLASSIFICATION.out.versions)
    }
    else if ( params.mode == 'donor_match' ){

        checkParams('demultiplexing_result', 'DONOR_MATCH', 'donor_match', true)

        ch_donor_match = ch_donor_match.map{
            meta, barcodes ->
            [meta, barcodes, params.demultiplexing_result]
        }

        if ( params.find_variants ){

            ['cell_genotype', 'vireo_filtered_variants'].each { p ->
                checkParams(p, 'FIND_VARIANTS', 'donor_match', true)
            }

            ch_find_variants = ch_find_variants.map{ meta ->
                [meta[0], params.cell_genotype, params.vireo_filtered_variants]
            }
        }

    }

    if (params.mode == 'genetic' | params.mode == 'hashing' | params.mode == 'rescue'){
        CREATE_ANNDATA_MUDATA(ch_create_anndata_mudata)
    }

    if (params.match_donor) {

        DONOR_MATCH(ch_donor_match,
            params.match_donor_method1 ?: [],
            params.match_donor_method2 ?: []
        )

        // there only is a best_intersect_assignment_after_match output in donor_match and rescue mode to run FIND_VARIANTS
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

            // subset gt_donors vcf with representative_variants
            // only vireo can produce gt_donors in rescue mode or user has to provide gt_donors in donor_match mode
            if (
                (params.mode == 'rescue' && params.genetic_tools && params.genetic_tools.split(',').contains('vireo')) |
                (params.mode == 'donor_match' && params.gt_donors && checkParams('gt_donors', 'SUBSET_GT_DONORS', 'donor_match', true))
            ) {

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

                ch_subset_gt_donors = ch_subset_gt_donors
                    .combine(params.mode == 'rescue'
                                ? GENETIC_DEMULTIPLEXING.out.gt_donors
                                : ch_subset_gt_donors.first().map{ meta, _variants, _type -> [meta, params.gt_donors] }
                                , by: 0)
                    .combine(DONOR_MATCH.out.best_donor_match, by: 0)

                SUBSET_GT_DONORS(ch_subset_gt_donors)

                ch_versions = ch_versions.mix(SUBSET_GT_DONORS.out.versions)
            }
            ch_versions = ch_versions.mix(FIND_VARIANTS.out.versions)
        }
        ch_versions = ch_versions.mix(DONOR_MATCH.out.versions)
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
