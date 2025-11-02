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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HADGE {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

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

    // TODO remove completely
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

    // TODO check again if they are the correct barcodes
    ch_donor_match = ch_genetic.map { meta, _rna, _hto, _bam, barcodes, _vcf ->
        [meta, barcodes]
    }

    // ------------------------------- preprocessing end --------------------------------


    //TODO mayybe go back to the orignal style with if, if else and so on
    //TODO geht hash und hash mode mit cell.genozyp?
    if (params.mode == 'genetic' || params.mode == 'rescue') {
        GENETIC_DEMULTIPLEXING(
            ch_genetic,
            params.genetic_tools.split(','),
            params.bam_qc,
            params.common_variants
        )

        ch_donor_match = ch_donor_match
            .join(GENETIC_DEMULTIPLEXING.out.summary_assignment)
            .join(GENETIC_DEMULTIPLEXING.out.cell_genotype)

        ch_versions = ch_versions.mix(GENETIC_DEMULTIPLEXING.out.versions)
    }

    if (params.mode == 'hashing' || params.mode == 'rescue') {
        HASH_DEMULTIPLEXING(
            ch_hashing,
            params.hash_tools.split(',')
        )

        ch_donor_match = ch_donor_match
            .join(HASH_DEMULTIPLEXING.out.summary_assignment)

        ch_versions = ch_versions.mix(HASH_DEMULTIPLEXING.out.versions)

        if(params.mode == 'rescue'){
            //TODO ext args for outer

            JOIN_RESULTS(ch_donor_match.map{
                meta, barcodes, gene_summary, cell_genotype, hash_summary ->
                [meta, [gene_summary,hash_summary]]
            })

            ch_donor_match = ch_donor_match
                .join(JOIN_RESULTS.out.csv)
                .map{
                    meta, barcodes, _gene_summary, cell_genotype, _hash_summary, joined_summary ->
                    [meta, barcodes, joined_summary, cell_genotype]
                }

            ch_versions = ch_versions.mix(JOIN_RESULTS.out.versions)
        }else{
            ch_donor_match = ch_donor_match
                .map{ meta, barcodes, hash_summary-> [meta, barcodes, hash_summary,[]] }
        }
    }

    if (params.mode == 'donor_match' || params.match_donor) {
        //TODO testing
        if (params.mode == 'donor_match'){
            ch_donor_match = ch_donor_match.map{
                meta, barcodes ->
                [meta, barcodes, params.demultiplexing_result, params.celldata, params.vireo_parent_dir]
            }
        }else{
            ch_donor_match.view()
            ch_donor_match = ch_donor_match.map{
                meta, barcodes, assignment_result, cell_genotype ->
                [meta, barcodes, assignment_result, cell_genotype, []]
            }
        }


        // TODO add params to nextflow.config and write DONOR matching
        DONOR_MATCH(ch_donor_match,
            params.match_donor_method1 ?: [],
            params.match_donor_method2 ?: [],
            params.findVariants,
            params.variant_count,
            params.variant_pct
        )


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
