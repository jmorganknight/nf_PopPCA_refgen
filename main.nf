nextflow.enable.dsl = 2

params.publish_dir_mode = params.publish_dir_mode ?: 'copy'

include { PARSE_MANIFEST }      from './modules/local/PARSE_MANIFEST'
include { NORMALIZE_INPUTS }    from './modules/local/NORMALIZE_INPUTS'
include { INTERSECT_VARIANTS }  from './modules/local/INTERSECT_VARIANTS'
include { QC_AND_PRUNE }        from './modules/local/QC_AND_PRUNE'
include { GENERATE_GLOBAL_PCA } from './modules/local/GENERATE_GLOBAL_PCA'
include { GENERATE_SUBLAYER_PCA } from './modules/local/GENERATE_SUBLAYER_PCA'
include { EMIT_PROVENANCE } from './modules/local/EMIT_PROVENANCE'

workflow {
    ref_manifest_file = file(params.ref_manifest, checkIfExists: true)
    ch_ref_manifest_file = channel.value(ref_manifest_file)
    ch_ref_manifest = channel.value(ref_manifest_file.toString())
    ch_params_yaml = channel.value(file("${projectDir}/assets/params.yaml", checkIfExists: true))
    ch_high_ld = channel.value(file(params.high_ld_regions, checkIfExists: true))

    def resolveRefPath = { String pathStr ->
        def p = new File(pathStr)
        if (p.isAbsolute()) {
            return file(pathStr, checkIfExists: true)
        }
        return file("${params.ref_dir}/${pathStr}", checkIfExists: true)
    }

    PARSE_MANIFEST(ch_ref_manifest)

    ch_manifest_record_list = PARSE_MANIFEST.out.active_datasets
        .splitText()
        .filter { line -> line?.trim() }
        .map { line -> new groovy.json.JsonSlurper().parseText(line) as Map }
        .collect()

    ch_active_datasets = ch_manifest_record_list
        .flatMap { records ->
            records.collect { meta ->
                meta.sample_count = (meta.sample_count ?: 0) as Integer
                tuple(meta, meta.paths.collect { pathItem -> resolveRefPath(pathItem.toString()) })
            }
        }

    ch_provenance_inputs = ch_manifest_record_list
        .map { records ->
            records
                .collectMany { meta -> meta.paths.collect { pathItem -> resolveRefPath(pathItem.toString()) } }
                .unique { staged -> staged.toString() }
        }

    NORMALIZE_INPUTS(ch_active_datasets)

    ch_intersection_inputs = NORMALIZE_INPUTS.out.normalized_panels
        .map { _meta, _pgen, _pvar, _psam, variantTable, _sampleTable, _summary -> variantTable }
        .collect()

    INTERSECT_VARIANTS(ch_intersection_inputs)

    ch_shared_variants = INTERSECT_VARIANTS.out.shared_variants

    QC_AND_PRUNE(
        NORMALIZE_INPUTS.out.normalized_panels,
        ch_shared_variants,
        ch_high_ld
    )

    ch_qc_records = QC_AND_PRUNE.out.qc_panels
        .map { meta, pgen, pvar, psam, sampleTable, summary ->
            [meta: meta, pgen: pgen, pvar: pvar, psam: psam, sample_table: sampleTable, summary: summary]
        }

    ch_global_inputs = ch_qc_records
        .collect()
        .map { records ->
            tuple(
                records.collect { record -> record.pgen },
                records.collect { record -> record.pvar },
                records.collect { record -> record.psam },
                records.collect { record -> record.sample_table },
                records.collect { record -> record.summary }
            )
        }

    GENERATE_GLOBAL_PCA(ch_global_inputs)

    ch_subgroup_inputs = ch_qc_records
        .collect()
        .map { records ->
            tuple(
                records.collect { record -> record.pgen },
                records.collect { record -> record.pvar },
                records.collect { record -> record.psam },
                records.collect { record -> record.sample_table },
                records.collect { record -> record.summary }
            )
        }

    GENERATE_SUBLAYER_PCA(ch_subgroup_inputs)

    EMIT_PROVENANCE(
        ch_ref_manifest_file,
        ch_provenance_inputs,
        ch_params_yaml,
        ch_high_ld,
        GENERATE_GLOBAL_PCA.out.summary,
        GENERATE_SUBLAYER_PCA.out.summary,
        GENERATE_GLOBAL_PCA.out.normalized_weights,
        channel.value((workflow.commitId ?: workflow.revision ?: 'unknown').toString()),
        channel.value(nextflow.version.toString()),
        channel.value((workflow.profile ?: '').toString()),
        channel.value((workflow.start ?: new Date()).toString()),
        channel.value((workflow.commandLine ?: '').toString()),
        channel.value((params.container ?: '').toString())
    )

    ch_final_manifest = GENERATE_GLOBAL_PCA.out.summary
        .mix(GENERATE_SUBLAYER_PCA.out.summary)
        .collect()
        .map { summaryFiles ->
            def parsed = summaryFiles.collect { summaryFile -> new groovy.json.JsonSlurper().parse(new File(summaryFile.toString())) as Map }
            def global = parsed.find { Map entry -> entry.layer == 'global' } ?: [:]
            def sublayerEntries = []
            parsed.findAll { Map entry -> entry.layer == 'sub' }.each { Map entry ->
                if (entry.models instanceof List) {
                    sublayerEntries.addAll(entry.models.collect { it as Map })
                } else {
                    sublayerEntries << entry
                }
            }
            def sublayers = sublayerEntries.sort { Map entry -> entry.superpop?.toString() ?: '' }
            def checksums = [:]
            parsed.each { entry ->
                (entry.files ?: [:]).each { key, value ->
                    checksums[key.toString()] = value.toString()
                }
                if (entry.models instanceof List) {
                    entry.models.each { modelEntry ->
                        (modelEntry.files ?: [:]).each { key, value ->
                            checksums[key.toString()] = value.toString()
                        }
                    }
                }
            }

            def outputManifest = [
                pipeline           : 'nf_PopPCA_refgen',
                target_build       : params.target_build,
                total_sample_count : (global.sample_count ?: 0) as Integer,
                generated_layers   : [
                    global   : [sample_count: (global.sample_count ?: 0) as Integer, variant_count: (global.variant_count ?: 0) as Integer],
                    sublayers: sublayers.collect { Map entry -> [superpop: entry.superpop, sample_count: (entry.sample_count ?: 0) as Integer, variant_count: (entry.variant_count ?: 0) as Integer] }
                ],
                file_checksums     : checksums
            ]
            groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(outputManifest))
        }

    ch_final_manifest.collectFile(name: 'reference_manifest.json', storeDir: "${params.outdir}/models", newLine: false)
}