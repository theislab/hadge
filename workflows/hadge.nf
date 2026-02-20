/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                    } from '../subworkflows/local/utils_nfcore_hadge_pipeline'
include { UNTAR as UNTAR_RNA                        } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_HTO                        } from '../modules/nf-core/untar/main'
include { EXTRACT_HASHES                            } from '../modules/local/extract_hashes/main'
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

workflow HADGE {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fasta          // file: /path/to/genome.fasta

    main:

    ch_versions = channel.empty()
    ch_donor_match = channel.empty()
    ch_find_variants = channel.empty()
    ch_subset_gt_donors = channel.empty()
    ch_multiqc_files = channel.empty()

    // ------------------------------ preprocessing start -------------------------------
    // untar matrices
    ch_rna = ch_samplesheet.map { meta, rna, _hto, _bam, _barcodes, _vcf -> [meta, rna] }
                    .branch { _meta, rna ->
                        tar: rna != null && rna.endsWith('.tar.gz')
                        directory: true
                    }


    ch_hto = ch_samplesheet.map { meta, _rna, hto, _bam, _barcodes, _vcf -> [meta, hto] }
                    .branch { _meta, hto ->
                        tar: hto != null && hto.endsWith('.tar.gz')
                        directory: true
                    }

    UNTAR_RNA(ch_rna.tar)
    ch_versions = ch_versions.mix(UNTAR_RNA.out.versions)

    UNTAR_HTO(ch_hto.tar)
    ch_versions = ch_versions.mix(UNTAR_HTO.out.versions)

    ch_rna = ch_rna.directory.mix(UNTAR_RNA.out.untar)
    ch_hto = ch_hto.directory.mix(UNTAR_HTO.out.untar)

    // extract hto names (hto can be null in genetic or donor_match mode)
    ch_hashes_non_null = EXTRACT_HASHES(ch_hto.filter { _meta, hto -> hto != null })
    ch_hashes_null = ch_hto.filter { _meta, hto -> hto == null }
    ch_hashes = ch_hashes_non_null.mix(ch_hashes_null)

    // join preprocessed channels
    ch_remaining_input = ch_samplesheet.map { meta, _rna, _hto, bam, barcodes, vcf -> [meta, bam, barcodes, vcf] }
    ch_preprocessed = ch_samplesheet.map { meta, _rna, _hto, _bam, _barcodes, _vcf -> [meta] }
                        .join(ch_rna)
                        .join(ch_hto)
                        .join(ch_remaining_input)
                        .join(ch_hashes)
                        .map {meta, rna, hto, bam, barcodes, vcf, hashes ->
                        if(hashes!= null){ meta += [hto_names: file(hashes).text.trim()] }
                        [meta, rna, hto, bam, barcodes, vcf]
                        }

    // create channels for deconvolution tools
    ch_genetic = ch_preprocessed.map { meta, rna, _hto, bam, barcodes, vcf -> [meta, bam, barcodes, vcf] }
    ch_hashing = ch_preprocessed.map { meta, rna, hto, _bam, _barcodes, _vcf -> [meta, rna, hto] }
    ch_create_anndata_mudata = ch_preprocessed.map { meta, rna, hto, _bam, _barcodes, _vcf -> [meta, rna, hto] }

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

        ch_donor_match = GENETIC_DEMULTIPLEXING.out.summary_assignment

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

        ch_donor_match = HASH_DEMULTIPLEXING.out.summary_assignment

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

        ch_donor_match = JOIN_RESULTS_ASSIGNMENT.out.csv

        if ( params.find_variants ){
            ch_find_variants = GENETIC_DEMULTIPLEXING.out.gt_cells
                .join(GENETIC_DEMULTIPLEXING.out.vireo_filtered_variants)
        }

        ch_versions = ch_versions.mix(GENETIC_DEMULTIPLEXING.out.versions)
        ch_versions = ch_versions.mix(HASH_DEMULTIPLEXING.out.versions)
        ch_versions = ch_versions.mix(JOIN_RESULTS_ASSIGNMENT.out.versions)
        ch_versions = ch_versions.mix(JOIN_RESULTS_CLASSIFICATION.out.versions)
    }
    else if ( params.mode == 'donor_match' ){

        ch_donor_match = ch_preprocessed.map{ meta, _rna, _hto, _bam, _barcodes, _vcf -> [meta, params.demultiplexing_result] }

        if ( params.find_variants ){
            ch_find_variants = ch_preprocessed.map{ meta, _rna, _hto, _bam, _barcodes, _vcf ->
                [meta, params.cell_genotype, params.vireo_filtered_variants]
            }
        }
    }

    if (params.mode == 'genetic' | params.mode == 'hashing' | params.mode == 'rescue'){
        CREATE_ANNDATA_MUDATA(
            ch_create_anndata_mudata.map { tuple ->
                // hto can be null in genetic mode
                if (params.mode == 'genetic'){ tuple.collect { it == null ? [] : it } }
                else{ tuple }
            }
        )
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
                .join(ch_donor_match)

            FIND_VARIANTS(
                ch_find_variants,
                params.variant_count,
                params.variant_pct
            )

            // subset gt_donors vcf with representative_variants
            if ( params.subset_gt_donors ) {
                ch_subset_gt_donors = FIND_VARIANTS.out.donor_specific_variants
                    .map { meta, subset_variants ->
                        tuple(meta, subset_variants, 'donor_specific')
                    }
                    .mix(
                        FIND_VARIANTS.out.vireo_variants
                            .map { meta, subset_variants ->
                                tuple(meta, subset_variants, 'vireo')
                            }
                    )

                ch_subset_gt_donors = params.mode == 'rescue'
                    ? ch_subset_gt_donors
                        .combine(GENETIC_DEMULTIPLEXING.out.gt_donors, by: 0)
                    : ch_subset_gt_donors
                        .map { meta, variants, type ->
                            [ meta, variants, type, params.gt_donors ]
                        }

                ch_subset_gt_donors = ch_subset_gt_donors
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
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
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
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

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
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}
