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

    main:
    ch_versions = Channel.empty()
    ch_vireo = Channel.empty()
    ch_demuxlet = Channel.empty()
    ch_freemuxlet = Channel.empty()
    ch_souporcell = Channel.empty()
    ch_cellsnp = Channel.empty()

    ch_summary = ch_samplesheet.map{ meta, rna, hto, _bam, barcodes, _vcf ->
        [meta, rna, hto, barcodes]
    }

    ch_samplesheet = ch_samplesheet.map{ meta, _rna, _hto, bam, barcodes, vcf ->
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


    if (methods.contains('vireo')) {
        SAMTOOLS_INDEX(ch_samplesheet.map { meta, bam, _barcodes, _vcf -> [meta, bam] })
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        CELLSNP_MODEA(
            ch_samplesheet.join(SAMTOOLS_INDEX.out.bai).map { meta, bam, barcodes, vcf, bai -> [meta, bam, bai, vcf, barcodes] }
        )

        ch_cellsnp = ch_cellsnp.mix(CELLSNP_MODEA.out.cell)
        ch_versions = ch_versions.mix(CELLSNP_MODEA.out.versions)

        VIREO(
            ch_samplesheet.join(CELLSNP_MODEA.out.cell).map { meta, _bam, _barcodes, vcf, cell -> [meta, cell, meta.n_samples, vcf, []] }
        )
        ch_vireo = ch_vireo.mix(VIREO.out.donor_ids)
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

        ch_soup_fasta = ch_samplesheet.map { meta, _bam, _barcodes, _vcf ->
            [ meta, file(params.ref) ]
        }

        // ch_soup_fasta = ch_samplesheet.map { meta, _bam, _barcodes, _vcf ->
        //     [ meta, file(params.fasta) ]
        // }

        //TODO update souporcell so that the first inputs also have the number of clusters

        SOUPORCELL(
            ch_souporcell_bam_barcodes_clusters,
            ch_soup_fasta
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

    ch_summary.view()

    GENE_SUMMARY(
        ch_summary,
        tuple(params.generate_anndata, params.generate_mudata)
    )

    ch_versions = ch_versions.mix(GENE_SUMMARY.out.versions)


    emit:
    summary_assignment = GENE_SUMMARY.out.assignment
    summary_classification = GENE_SUMMARY.out.classification
    cell_genotype = ch_cellsnp
    versions = ch_versions // channel: [ versions.yml ]
}
