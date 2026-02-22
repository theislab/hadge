include { BAM_QC             } from '../bam_qc'
include { FILTER_BAM         } from '../../../modules/local/filter_bam'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index'
include { CELLSNP_MODEA      } from '../../../modules/nf-core/cellsnp/modea'
include { VIREO              } from '../../../modules/nf-core/vireo'
include { POPSCLE_DSCPILEUP  } from '../../../modules/nf-core/popscle/dscpileup'
include { POPSCLE_DEMUXLET   } from '../../../modules/nf-core/popscle/demuxlet'
include { POPSCLE_FREEMUXLET } from '../../../modules/nf-core/popscle/freemuxlet'
include { SOUPORCELL         } from '../../../modules/nf-core/souporcell'
include { GENE_SUMMARY       } from '../../../modules/local/gene_summary'

workflow GENETIC_DEMULTIPLEXING {
    take:
    ch_samplesheet  // channel: samplesheet read in from --input
    methods         // list of strings
    bam_qc          // boolean
    common_variants // file
    fasta           // file: /path/to/genome.fasta

    main:
    ch_versions = channel.empty()
    ch_vireo = channel.empty()
    ch_demuxlet = channel.empty()
    ch_freemuxlet = channel.empty()
    ch_souporcell = channel.empty()
    ch_gt_cells = channel.empty()
    ch_gt_donors =  channel.empty()
    ch_vireo_filtered_variants = channel.empty()

    ch_summary = ch_samplesheet.map{ meta, _bam, barcodes, _vcf ->
        [meta, barcodes]
    }

    ch_samplesheet = ch_samplesheet.map{ meta, bam, barcodes, vcf ->
        [meta, bam, barcodes, vcf]
    }

    if (bam_qc) {
        BAM_QC(ch_samplesheet.map { meta, bam, _barcodes, _vcf -> [meta, bam] })
        ch_versions = ch_versions.mix(BAM_QC.out.versions)

        ch_samplesheet = ch_samplesheet
            .join(BAM_QC.out.bam)
            .map { meta, _bam, barcodes, vcf, new_bam -> [meta, new_bam, barcodes, vcf] }
    }

    if (common_variants) {
        FILTER_BAM(
            ch_samplesheet.map { meta, bam, barcodes, _vcf ->
                [
                    meta,
                    bam,
                    barcodes,
                ]
            },
            common_variants,
        )
        ch_versions = ch_versions.mix(FILTER_BAM.out.versions)

        ch_samplesheet = ch_samplesheet
            .join(FILTER_BAM.out.bam)
            .map { meta, _bam, barcodes, vcf, new_bam -> [meta, new_bam, barcodes, vcf] }
    }



    if ( params.find_variants | methods.contains('vireo')){
        SAMTOOLS_INDEX(ch_samplesheet.map { meta, bam, _barcodes, _vcf -> [meta, bam] })
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        CELLSNP_MODEA(
            ch_samplesheet.join(SAMTOOLS_INDEX.out.bai).map { meta, bam, barcodes, vcf, bai -> [meta, bam, bai, vcf, barcodes] }
        )

        ch_gt_cells = ch_gt_cells.mix(CELLSNP_MODEA.out.cell)
        ch_versions = ch_versions.mix(CELLSNP_MODEA.out.versions)
    }

    if (methods.contains('vireo')) {
        VIREO(
            ch_samplesheet.join(ch_gt_cells).map { meta, _bam, _barcodes, vcf, cell -> [meta, cell, meta.n_samples, vcf, []] }
        )

        ch_vireo = ch_vireo.mix(VIREO.out.donor_ids)
        ch_vireo_filtered_variants = ch_vireo_filtered_variants.mix(VIREO.out.filtered_variants)
        ch_gt_donors = ch_gt_donors.mix(VIREO.out.genotype_vcf)
        ch_versions = ch_versions.mix(VIREO.out.versions)
    }

    if (methods.contains('demuxlet') || methods.contains('freemuxlet')) {
        ch_dscpileup = ch_samplesheet.map { meta, bam, _barcodes, vcf -> [meta, bam, vcf] }
        POPSCLE_DSCPILEUP(ch_dscpileup)
        ch_versions = ch_versions.mix(POPSCLE_DSCPILEUP.out.versions)

        if (methods.contains('demuxlet')) {
            ch_demuxlet_input = POPSCLE_DSCPILEUP.out.plp.join(ch_samplesheet).map { meta, plp, bam, _barcodes, vcf -> [meta, plp, bam, vcf] }
            POPSCLE_DEMUXLET(ch_demuxlet_input)
            ch_demuxlet = ch_demuxlet.mix(POPSCLE_DEMUXLET.out.demuxlet_result)
            ch_versions = ch_versions.mix(POPSCLE_DEMUXLET.out.versions)

        }
        if (methods.contains('freemuxlet')) {
            ch_freemuxlet_input = POPSCLE_DSCPILEUP.out.directory.join(ch_samplesheet).map { meta, plp_dir, _bam, _barcodes, _vcf -> [meta, plp_dir, meta.n_samples] }
            POPSCLE_FREEMUXLET(ch_freemuxlet_input)
            ch_freemuxlet = ch_freemuxlet.mix(POPSCLE_FREEMUXLET.out.result)
            ch_versions = ch_versions.mix(POPSCLE_FREEMUXLET.out.versions)
        }
    }

    if (methods.contains('souporcell')) {

        ch_souporcell_bam_barcodes_clusters = ch_samplesheet.map { meta, bam, barcodes, _vcf ->
            [ meta, bam, barcodes, meta.n_samples ]
        }

        if (! params.fasta ) {
            log.warn("The pipeline is downloading the full reference genome from ${params.genome} because only `genome` and not `fasta` is set. To reduce long download times and high bandwidth usage, provide your own reference by specifying `fasta`.")
        }

        SOUPORCELL(
            ch_souporcell_bam_barcodes_clusters,
            channel.value([[id: 'fasta'], file(fasta, checkIfExists: true)])
        )

        ch_souporcell = ch_souporcell.mix(SOUPORCELL.out.clusters)
        ch_versions = ch_versions.mix(SOUPORCELL.out.versions)
    }

    ch_summary = ch_summary
        .join(ch_vireo, remainder: true)
        .join(ch_demuxlet, remainder: true)
        .join(ch_freemuxlet, remainder: true)
        .join(ch_souporcell, remainder: true)
        .map { tuple -> tuple.collect { it == null ? [] : it } }

    GENE_SUMMARY(ch_summary)

    ch_versions = ch_versions.mix(GENE_SUMMARY.out.versions)


    emit:
    summary_assignment = GENE_SUMMARY.out.assignment
    summary_classification = GENE_SUMMARY.out.classification
    vireo_filtered_variants = ch_vireo_filtered_variants
    gt_cells = ch_gt_cells
    gt_donors = ch_gt_donors
    versions = ch_versions // channel: [ versions.yml ]
}
