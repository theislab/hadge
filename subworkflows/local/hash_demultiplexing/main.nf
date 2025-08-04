include { UNTAR as UNTAR_RNA                                       } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_HTO                                       } from '../../../modules/nf-core/untar'
include { RENAME_GENES_TO_FEATURES as RENAME_GENES_TO_FEATURES_RNA } from '../../../modules/local/rename_genes_to_features'
include { RENAME_GENES_TO_FEATURES as RENAME_GENES_TO_FEATURES_HTO } from '../../../modules/local/rename_genes_to_features'
include { DROPLETUTILS_MTXCONVERT as MTXCONVERT_RNA                } from '../../../modules/local/dropletutils/mtxconvert'
include { DROPLETUTILS_MTXCONVERT as MTXCONVERT_HTO                } from '../../../modules/local/dropletutils/mtxconvert'
include { PREPROCESSING_FOR_HTODEMUX_MULTISEQ                      } from '../../../modules/local/preprocessing_for_htodemux_multiseq'
include { HTODEMUX                                                 } from '../../../modules/nf-core/htodemux'
include { HTODEMUX_VISUALIZATION                                   } from '../../../modules/local/htodemux_visualization'
include { MULTISEQDEMUX                                            } from '../../../modules/nf-core/multiseqdemux'
include { DEMUXEM                                                  } from '../../../modules/nf-core/demuxem'
include { GMMDEMUX                                                 } from '../../../modules/nf-core/gmmdemux'
include { HASHEDDROPS                                              } from '../../../modules/nf-core/hasheddrops'
include { HASH_SUMMARY                                             } from '../../../modules/local/hash_summary'

workflow HASH_DEMULTIPLEXING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    methods        // list of strings

    main:

    ch_results = Channel.empty()

    ch_versions = Channel.empty()

    ch_htodemux_assignments = Channel.empty()
    ch_htodemux_classifications = Channel.empty()
    ch_multiseq = Channel.empty()
    ch_cellhashr = Channel.empty()
    ch_demuxem = Channel.empty()
    ch_gmmdemux = Channel.empty()
    ch_hasheddrops = Channel.empty()
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

    ch_rna = ch_samplesheet.map { meta, rna, _hto -> [meta, rna] }
                    .branch { _meta, rna ->
                        tar: rna.endsWith('.tar.gz')
                        directory: true
                    }
    ch_hto = ch_samplesheet.map { meta, _rna, hto -> [meta, hto] }
                    .branch { _meta, hto ->
                        tar: hto.endsWith('.tar.gz')
                        directory: true
                    }

    UNTAR_RNA(ch_rna.tar)
    ch_versions = ch_versions.mix(UNTAR_RNA.out.versions)

    UNTAR_HTO(ch_hto.tar)
    ch_versions = ch_versions.mix(UNTAR_HTO.out.versions)

    ch_rna = ch_rna.directory.mix(UNTAR_RNA.out.untar)
    ch_hto = ch_hto.directory.mix(UNTAR_HTO.out.untar)

    ch_rna = RENAME_GENES_TO_FEATURES_RNA(ch_rna)
    ch_hto = RENAME_GENES_TO_FEATURES_HTO(ch_hto)

    ch_samplesheet = ch_samplesheet.map { meta, _rna, _hto -> [meta] }.join(ch_rna).join(ch_hto)

    if (methods.contains('htodemux') || methods.contains('multiseq')) {
        PREPROCESSING_FOR_HTODEMUX_MULTISEQ(
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(PREPROCESSING_FOR_HTODEMUX_MULTISEQ.out.versions)

        if (methods.contains('htodemux')) {
            HTODEMUX(
                PREPROCESSING_FOR_HTODEMUX_MULTISEQ.out.seurat_object.map { meta, seurat_object -> [meta, seurat_object, "HTO"] }
            )
            def out = HTODEMUX.out

            //HTODEMUX.out.assignment.view()
            // ch_results = ch_results.mix(
            //     HTODEMUX.out.assignment
            //     .join(HTODEMUX.out.classification)
            //     .join(HTODEMUX.out.params)
            //     .map { meta, assign, classi, par ->
            //         def files = [assignment: assign, classification: classi, params: par]
            //         [meta, [results: files, method: 'htodemux']]
            //     })

            // [meta, [result: assign, method:'htodemux_assignment' ]
            // [meta, [result: classi, method:'htodemux_assignment' ]


            // assignmethod:, 'htodemux_classification'


            // ch_versions = ch_versions.mix(HTODEMUX.out.versions)

            ch_assignments = HTODEMUX.out.assignment
                .map { meta, assignment ->
                    [meta, [result: assignment, method: 'htodemux_assignment']]
                }

            ch_classifications = HTODEMUX.out.classification
                .map { meta, classification ->
                    [meta, [result: classification, method: 'htodemux_classification']]
                }

            ch_results = ch_results
                .mix(ch_assignments,ch_classifications)

            println("results oben ---->")
            //ch_results.view()
            println("results oben ---->")


            //ch_results = ch_results.mix(HTODEMUX.out.assignment.map { meta, result -> [meta, [path: result, method: 'htodemux']] })

            //HTODEMUX.out.view()
            //HTODEMUX.out.assignment.join(HTODEMUX.out.classification).join(HTODEMUX.out.params).view("testiiii "+ it)


            ch_htodemux_assignments = ch_htodemux_assignments.mix(HTODEMUX.out.assignment)
            ch_htodemux_classifications = ch_htodemux_classifications.mix(HTODEMUX.out.classification)



            HTODEMUX_VISUALIZATION(
                HTODEMUX.out.rds.map { meta, seurat_object -> [meta, seurat_object, "HTO"] }
            )
            ch_versions = ch_versions.mix(HTODEMUX_VISUALIZATION.out.versions)
        }
        if (methods.contains('multiseq')) {
            MULTISEQDEMUX(
                PREPROCESSING_FOR_HTODEMUX_MULTISEQ.out.seurat_object.map { meta, seurat_object -> [meta, seurat_object, "HTO"] }
            )
            // MULTISEQDEMUX.out.results.view()
            // ch_results = ch_results.mix(
            //     MULTISEQDEMUX.out.results
            //     .map { meta, results ->
            //         def paths = [assignment: results, classification: null]
            //         [meta, [results: paths, method: 'htodemux']]
            //     })

            ch_results = ch_results.mix(
                MULTISEQDEMUX.out.results
                .map { meta, file ->
                    [meta, [results: file, method: 'multiseq']]
                })


            //ch_results = ch_results.mix(MULTISEQDEMUX.out.results.map { meta, result -> [meta, [path: result, method: 'multiseq']] })
            ch_multiseq = ch_multiseq.mix(MULTISEQDEMUX.out.results)
            ch_versions = ch_versions.mix(MULTISEQDEMUX.out.versions)
        }
    }

    // TODO rename to bff since we named the module bff
    if (methods.contains('cellhashr')) {
        error("CellHashR not implemented")
    }
    if (methods.contains('demuxem')) {
        ch_samplesheet.map { meta, rna, hto ->
            {
                if (!rna) {
                    error("RNA matrix not provided for sample ${meta.id}, but this is required for DemuxEM. Please check your input samplesheet.")
                }
                if (!hto) {
                    error("HTO matrix not provided for sample ${meta.id}, but this is required for DemuxEM. Please check your input samplesheet.")
                }
            }
        }

        MTXCONVERT_RNA(ch_samplesheet.map { meta, rna, _hto -> [meta, rna] }, false)
        ch_versions = ch_versions.mix(MTXCONVERT_RNA.out.versions)

        MTXCONVERT_HTO(ch_samplesheet.map { meta, _rna, hto -> [meta, hto] }, true)
        ch_versions = ch_versions.mix(MTXCONVERT_HTO.out.versions)

        DEMUXEM(
            MTXCONVERT_RNA.out.h5.join(MTXCONVERT_HTO.out.csv),
            params.demuxem_gender_genes,
            params.genome ?: [],
            true,
        )

        ch_demuxem = ch_demuxem.mix(DEMUXEM.out.out_zarr)
        ch_versions = ch_versions.mix(DEMUXEM.out.versions)
    }
    if (methods.contains('gmm-demux')) {
        ch_gmmdemux_input = ch_samplesheet.map { meta, _rna, hto -> [meta, hto, "MS-11,MS-12", meta.n_cells] }

        ch_gmmdemux_input.map { meta, hto, hto_names, _estimated_cells ->
            {
                if (!hto) {
                    error("HTO matrix not provided for sample ${meta.id}, but this is required for GMM-Demux. Please check your input samplesheet.")
                }
                if (!hto_names) {
                    error("HTO names not provided for sample ${meta.id}, but this is required for GMM-Demux. Please check your input samplesheet.")
                }
            }
        }
        GMMDEMUX(
            ch_gmmdemux_input,
            true,
            true,
            [],
            [],
        )
        ch_versions = ch_versions.mix(GMMDEMUX.out.versions)
    }
    if (methods.contains('hasheddrops')) {
        HASHEDDROPS(
            ch_samplesheet.map { meta, rna, hto -> [meta, hto, "FALSE", rna] }
        )
        ch_hasheddrops = ch_hasheddrops.mix(HASHEDDROPS.out.results)
        ch_versions = ch_versions.mix(HASHEDDROPS.out.versions)
    }
    if (methods.contains('hashsolo')) {
        error("HashSolo not implemented")
    }

    //ch_results.view()
    // group by module, give the channel to summary module
    // group by module
    //ch_results_grouped = ch_results.groupTuple(by: 0).view()



    def methods_list = methods as List
    def sorted_method_files = ['htodemux_assignment', 'htodemux_classification', 'multiseq', 'cellhashr', 'demuxem', 'gmm-demux', 'hasheddrops', 'hashsolo']

    def used_method_files = methods_list.contains('htodemux')
        ? (methods_list - 'htodemux') + ['htodemux_assignment', 'htodemux_classification']
        : methods_list

    def empty_method_files = sorted_method_files - used_method_files

    println("Empty" + empty_method_files)
    println("Used" + used_method_files)

    // ch_hashing_summary = ch_results.groupTuple(by: 0).view()
    //ch_results.view{"results unten"+it}
    // def met = ['htodemux', 'multiseq']
    // def diff_met = sorted_methods - met
    // def diff_methods = sorted_methods - used_methods
    // println("methods "+methods)
    // println("met "+met)
    // println("methods type "+ methods.getClass())
    // println("met "+met.getClass())
    // println("sorted - methods "+diff_methods)
    // println("sorted - met "+diff_met)
    // sort the methods result paths as in sorted_methods and add empty results for methods not calculated
    // you either have a single file path, a list (groovy map) of file paths or null if there where no results
    // e.g. [meta, file1, file2, [A: file3, B: file4], null, file5, ...]

    //ch_results_sorted = ch_results.groupTuple(by: 0).view()

    ch_results_sorted = ch_results.groupTuple(by: 0)
    .map { meta, results ->
        def empty_results = empty_method_files.collect {[results: null,method: it] }
        // println("empties: "+ empty_results)
        // println("+ --> "+(results + empty_results))
        def sorted_results = (results + empty_results)
            .sort { a, b ->
    sorted_method_files.indexOf(a.method) <=> sorted_method_files.indexOf(b.method)
 }
            .collect { it.results }
        // println("sorted: "+ sorted_results)
        [meta] + sorted_results
    }
    // .view{
    //     "final "+ it
    // }


    // ch_samplesheet.join(ch_results_sorted).view()


    def generate_anndata = false
    def generate_mudata = false

    // TODO

    ['htodemux_assignment', 'htodemux_classification', 'multiseq', 'cellhashr', 'demuxem', 'gmm-demux', 'hasheddrops', 'hashsolo']

    ch_summary = ch_samplesheet
        .join(ch_htodemux_assignments, remainder: true)
        .join(ch_htodemux_classifications, remainder: true)
        .join(ch_multiseq, remainder: true)
        .join(ch_cellhashr, remainder: true)
        .join(ch_demuxem , remainder: true)
        .join(ch_gmmdemux, remainder: true)
        .join(ch_hasheddrops, remainder: true)
        .join(ch_hashsolo, remainder: true)
        .map { tuple -> tuple.collect { it == null ? [] : it } }
    // Empty inputs solved as recommended here:
    // https://nf-co.re/docs/guidelines/components/modules#optional-inputs

    ch_summary.view()

    HASH_SUMMARY(ch_summary, generate_anndata, generate_mudata)



    emit:
    versions = ch_versions // channel: [ versions.yml ]
}
