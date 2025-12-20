include { PREPROCESSING_FOR_HTODEMUX_MULTISEQ                      } from '../../../modules/local/preprocessing_for_htodemux_multiseq'
include { HTODEMUX                                                 } from '../../../modules/nf-core/htodemux'
include { HTODEMUX_VISUALIZATION                                   } from '../../../modules/local/htodemux_visualization'
include { MULTISEQDEMUX                                            } from '../../../modules/nf-core/multiseqdemux'
include { BFF                                                      } from '../../../modules/nf-core/bff'
include { DEMUXEM                                                  } from '../../../modules/nf-core/demuxem'
include { DROPLETUTILS_MTXCONVERT as MTXCONVERT_RNA                } from '../../../modules/local/dropletutils/mtxconvert/main'
include { DROPLETUTILS_MTXCONVERT as MTXCONVERT_HTO                } from '../../../modules/local/dropletutils/mtxconvert/main'
include { GMMDEMUX                                                 } from '../../../modules/nf-core/gmmdemux'
include { SCANPY_HASHSOLO as HASHSOLO                              } from '../../../modules/nf-core/scanpy/hashsolo'
include { HASHEDDROPS                                              } from '../../../modules/nf-core/hasheddrops'
include { HASH_SUMMARY                                             } from '../../../modules/local/hash_summary'


workflow HASH_DEMULTIPLEXING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    methods        // list of strings

    main:

    ch_versions = Channel.empty()

    ch_htodemux_assignments = Channel.empty()
    ch_htodemux_classifications = Channel.empty()
    ch_multiseq = Channel.empty()
    ch_bff = Channel.empty()
    ch_demuxem = Channel.empty()
    ch_gmmdemux_results = Channel.empty()
    ch_gmmdemux_config = Channel.empty()
    ch_hasheddrops_results = Channel.empty()
    ch_hasheddrops_id_to_hash = Channel.empty()
    ch_hashsolo = Channel.empty()

    ch_samplesheet.map { meta, rna, hto ->
        {
            if (!rna) {
                error("RNA matrix not provided for sample ${meta.id}, but this is required for hash demultiplexing. Please check your input samplesheet.")
            }
            if (!hto) {
                error("HTO matrix not provided for sample ${meta.id}, but this is required for hash demultiplexing. Please check your input samplesheet.")
            }
        }
    }

    if (methods.contains('htodemux') || methods.contains('multiseq')) {
        ch_samplesheet.map { meta, rna, hto ->
            if(meta.hto_names.split(",").any { it.contains('_') }){
                def bad = meta.hto_names.split(",").findAll { it.contains('_') }.join(', ')
                throw new IllegalArgumentException(
                    "Running hadge with the methods htodemux or multiseq does not allow to use underscores ('_') in HTO names. Both tools require a SeuratObject as input, which will replace '_' with '-' leading to ambiguous or misleading assignment summaries. Please remove underscores ('_') from: ${bad}"
                )
            }
        }

        PREPROCESSING_FOR_HTODEMUX_MULTISEQ(
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(PREPROCESSING_FOR_HTODEMUX_MULTISEQ.out.versions)

        if (methods.contains('htodemux')) {
            HTODEMUX(
                PREPROCESSING_FOR_HTODEMUX_MULTISEQ.out.seurat_object.map { meta, seurat_object -> [meta, seurat_object, "HTO"] }
            )

            ch_assignments = HTODEMUX.out.assignment
                .map { meta, assignment ->
                    [meta, [result: assignment, method: 'htodemux_assignment']]
                }

            ch_classifications = HTODEMUX.out.classification
                .map { meta, classification ->
                    [meta, [result: classification, method: 'htodemux_classification']]
                }

            ch_htodemux_assignments = ch_htodemux_assignments.mix(HTODEMUX.out.assignment)
            ch_htodemux_classifications = ch_htodemux_classifications.mix(HTODEMUX.out.classification)

            HTODEMUX_VISUALIZATION(
                HTODEMUX.out.rds.map { meta, seurat_object -> [meta, seurat_object, "HTO"] }
            )
            ch_versions = ch_versions.mix(HTODEMUX_VISUALIZATION.out.versions)
        }
        if (methods.contains('multiseq')) {
            MULTISEQDEMUX(
                PREPROCESSING_FOR_HTODEMUX_MULTISEQ.out.seurat_object.map { meta, seurat_object -> [
                        meta,
                        seurat_object,
                        "HTO"
                    ]
                }
            )

            ch_multiseq = ch_multiseq.mix(MULTISEQDEMUX.out.results)
            ch_versions = ch_versions.mix(MULTISEQDEMUX.out.versions)
        }
    }

    if (methods.contains('bff')) {
        BFF(ch_samplesheet.map { meta, _rna, hto -> [meta,hto,params.bff_methods,params.bff_preprocessing]})
        ch_bff = ch_bff.mix(BFF.out.assignment)
        ch_versions = ch_versions.mix(BFF.out.versions)
    }

    if (methods.contains('demuxem')) {

        MTXCONVERT_RNA(ch_samplesheet.map { meta, rna, _hto -> [meta, rna] }, false)
        ch_versions = ch_versions.mix(MTXCONVERT_RNA.out.versions)

        MTXCONVERT_HTO(ch_samplesheet.map { meta, _rna, hto -> [meta, hto] }, true)
        ch_versions = ch_versions.mix(MTXCONVERT_HTO.out.versions)

        DEMUXEM(
            MTXCONVERT_RNA.out.h5.join(MTXCONVERT_HTO.out.csv),
            params.demuxem_gender_genes,
            params.genome ?: [],
            params.demuxem_generate_diagnostic_plots
        )

        ch_demuxem = ch_demuxem.mix(DEMUXEM.out.out_zarr)
        ch_versions = ch_versions.mix(DEMUXEM.out.versions)
    }

    if (methods.contains('gmm-demux')) {

        ch_gmmdemux_input = ch_samplesheet.map { meta, _rna, hto -> [
                    meta,
                    hto,
                    params.gmmdemux_hto_names ? params.gmmdemux_hto_names : meta.hto_names,
                    params.gmmdemux_estimated_n_cells ? gmmdemux_estimated_n_cells : [],
                ]
            }

        GMMDEMUX(
            ch_gmmdemux_input,
            params.gmmdemux_type_report,
            params.gmmdemux_summary_report,
            params.gmmdemux_skip ? params.gmmdemux_skip : [],
            params.gmmdemux_examine ? params.gmmdemux_examine : []
        )

        ch_versions = ch_versions.mix(GMMDEMUX.out.versions)
        ch_gmmdemux_results = ch_gmmdemux_results.mix(GMMDEMUX.out.classification_report)
        ch_gmmdemux_config = ch_gmmdemux_config.mix(GMMDEMUX.out.config_report)
    }
    if (methods.contains('hasheddrops')) {

        HASHEDDROPS(
            ch_samplesheet.map { meta, rna, hto -> [
                    meta,
                    hto,
                    params.hasheddrops_runEmptyDrops.toString().toUpperCase(),
                    rna
                ]
            }
        )

        ch_hasheddrops_results = ch_hasheddrops_results.mix(HASHEDDROPS.out.results)
        ch_hasheddrops_id_to_hash = ch_hasheddrops_id_to_hash.mix(HASHEDDROPS.out.id_to_hash)
        ch_versions = ch_versions.mix(HASHEDDROPS.out.versions)
    }
    if (methods.contains('hashsolo')) {

        HASHSOLO(
            ch_samplesheet.map {meta, _rna, hto -> [
                    meta,
                    hto,
                    params.hashsolo_cell_hashing_columns ? params.hashsolo_cell_hashing_columns : []
                ]
            }
        )

        ch_hashsolo = ch_hashsolo.mix(HASHSOLO.out.assignment)
        ch_versions = ch_versions.mix(HASHSOLO.out.versions)
    }

    ch_summary = ch_samplesheet
        .join(ch_htodemux_assignments, remainder: true)
        .join(ch_htodemux_classifications, remainder: true)
        .join(ch_multiseq, remainder: true)
        .join(ch_bff, remainder: true)
        .join(ch_demuxem , remainder: true)
        .join(ch_gmmdemux_results, remainder: true)
        .join(ch_gmmdemux_config, remainder: true)
        .join(ch_hasheddrops_results, remainder: true)
        .join(ch_hasheddrops_id_to_hash, remainder: true)
        .join(ch_hashsolo, remainder: true)
        .map { tuple -> tuple.collect { it == null ? [] : it } }
    // Empty inputs solved as recommended here:
    // https://nf-co.re/docs/guidelines/components/modules#optional-inputs

    HASH_SUMMARY(ch_summary,params.bff_methods)

    ch_versions = ch_versions.mix(HASH_SUMMARY.out.versions)

    emit:
    summary_assignment = HASH_SUMMARY.out.assignment
    summary_classification = HASH_SUMMARY.out.classification
    versions = ch_versions // channel: [ versions.yml ]
}
