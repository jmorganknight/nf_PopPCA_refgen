process PARSE_MANIFEST {
    tag { "manifest:${new File(ref_manifest.toString()).name}" }

    input:
    val ref_manifest

    output:
    path 'active_datasets.jsonl', emit: active_datasets
    path 'manifest_summary.json', emit: summary

    exec:
    def manifestFile = new File(ref_manifest.toString())
    if (!manifestFile.exists()) {
        throw new IllegalArgumentException("Manifest path does not exist: ${ref_manifest}")
    }
    def payload = new groovy.yaml.YamlSlurper().parseText(manifestFile.text) as Map
    def datasets = (payload.datasets ?: []) as List<Map>
    def allowedAssays = ['WGS', 'WES', 'Microarray'] as Set
    def allowedFormats = ['vcf', 'plink1', 'plink2'] as Set
    def allowedSuperpops = ['GLOBAL', 'EUR', 'AFR', 'EAS', 'SAS', 'AMR', 'MENA'] as Set

    if (!datasets) {
        throw new IllegalArgumentException("No datasets were declared in ${ref_manifest}")
    }

    def parsePathEntry = { Object rawItem, String context ->
        if (rawItem instanceof Map) {
            def itemMap = rawItem as Map
            def parsedPath = (itemMap.path ?: '').toString().trim()
            if (!parsedPath) {
                throw new IllegalArgumentException("${context} is missing required key 'path'")
            }
            def expected = (itemMap.expected_sha256 ?: '').toString().trim()
            return [path: parsedPath, expected_sha256: expected]
        }

        def parsedPath = (rawItem ?: '').toString().trim()
        if (!parsedPath) {
            throw new IllegalArgumentException("${context} contains an empty path value")
        }
        return [path: parsedPath, expected_sha256: '']
    }

    def resolvePath = { String rawPath ->
        def p = new File(rawPath)
        if (p.isAbsolute()) {
            return p
        }
        return new File(params.ref_dir.toString(), rawPath)
    }

    def sha256File = { File f ->
        def proc = ["sha256sum", f.absolutePath].execute()
        def stdout = new StringBuffer()
        def stderr = new StringBuffer()
        proc.consumeProcessOutput(stdout, stderr)
        int rc = proc.waitFor()
        if (rc != 0) {
            throw new IllegalStateException("Failed computing sha256 for '${f.absolutePath}': ${stderr.toString().trim()}")
        }
        def output = stdout.toString().trim()
        if (!output) {
            throw new IllegalStateException("Failed computing sha256 for '${f.absolutePath}': empty output")
        }
        return output.tokenize()[0]
    }

    def seenIds = [] as Set
    def activeRecords = datasets.findAll { Map dataset -> (dataset.active ?: false) as Boolean }.collect { Map dataset ->
        def datasetId = (dataset.id ?: '').toString()
        if (!(datasetId ==~ /[A-Za-z0-9_.-]+/)) {
            throw new IllegalArgumentException("Dataset id '${datasetId}' must match [A-Za-z0-9_.-]+")
        }
        if (!seenIds.add(datasetId)) {
            throw new IllegalArgumentException("Duplicate dataset id '${datasetId}' detected in ${ref_manifest}")
        }

        def assay = (dataset.assay ?: '').toString()
        if (!allowedAssays.contains(assay)) {
            throw new IllegalArgumentException("Dataset '${datasetId}' has unsupported assay '${assay}'")
        }

        def superpop = (dataset.superpop ?: '').toString()
        if (!allowedSuperpops.contains(superpop)) {
            throw new IllegalArgumentException("Dataset '${datasetId}' has unsupported superpop '${superpop}'")
        }

        def build = (dataset.build ?: '').toString()
        if (!build) {
            throw new IllegalArgumentException("Dataset '${datasetId}' is missing build")
        }

        def input = (dataset.input ?: [:]) as Map
        def inputFormat = (input.format ?: '').toString()
        if (!allowedFormats.contains(inputFormat)) {
            throw new IllegalArgumentException("Dataset '${datasetId}' has unsupported input format '${inputFormat}'")
        }

        def pathEntries = ((input.paths ?: []) as List).withIndex().collect { Object pathItem, int i ->
            parsePathEntry(pathItem, "Dataset '${datasetId}' input.paths[${i}]")
        }
        if (!pathEntries) {
            throw new IllegalArgumentException("Dataset '${datasetId}' must provide one or more staged input paths")
        }

        def metadata = (dataset.metadata ?: [:]) as Map
        def metadataPath = ''
        def metadataExpectedSha256 = ''
        if (metadata.path != null && metadata.path.toString().trim()) {
            def parsedMetadataPath = parsePathEntry(metadata.path, "Dataset '${datasetId}' metadata.path")
            metadataPath = parsedMetadataPath.path
            metadataExpectedSha256 = parsedMetadataPath.expected_sha256
        }
        if (metadataPath) {
            pathEntries << [path: metadataPath, expected_sha256: metadataExpectedSha256]
        }

        pathEntries.each { Map pathEntry ->
            def expected = (pathEntry.expected_sha256 ?: '').toString().trim()
            if (!expected) {
                return
            }

            def inputFile = resolvePath(pathEntry.path.toString())
            if (!inputFile.exists() || !inputFile.isFile()) {
                throw new IllegalArgumentException("Dataset '${datasetId}' checksum verification failed: file not found for expected_sha256 at '${pathEntry.path}'")
            }

            def observed = sha256File(inputFile)
            if (!observed.equalsIgnoreCase(expected)) {
                throw new IllegalStateException("Dataset '${datasetId}' checksum mismatch for '${pathEntry.path}'. expected='${expected}' observed='${observed}'")
            }
        }

        def paths = pathEntries.collect { Map pathEntry -> pathEntry.path.toString() }

        def sampleIdCol = (metadata.sample_id_col ?: '#IID').toString()
        def popCol = (metadata.pop_col ?: '').toString()
        def superpopCol = (metadata.superpop_col ?: '').toString()

        [
            id          : datasetId,
            label       : (dataset.label ?: datasetId).toString(),
            assay       : assay,
            superpop    : superpop,
            sample_count: (dataset.sample_count ?: 0) as Integer,
            build       : build,
            input_format: inputFormat,
            paths       : paths,
            metadata_path: metadataPath,
            metadata_sample_id_col: sampleIdCol,
            metadata_pop_col: popCol,
            metadata_superpop_col: superpopCol,
        ]
    }

    if (!activeRecords) {
        throw new IllegalArgumentException("No active datasets were found in ${ref_manifest}")
    }

    def activeOut = new File(task.workDir.toString(), 'active_datasets.jsonl')
    activeOut.text = activeRecords.collect { Map record -> groovy.json.JsonOutput.toJson(record) }.join('\n') + '\n'

    def summary = [
        manifest_version     : payload.manifest_version ?: 1,
        target_build         : (payload.target_build ?: 'GRCh38').toString(),
        active_dataset_count : activeRecords.size(),
        active_sample_count  : activeRecords.collect { Map record -> record.sample_count as Integer }.sum() ?: 0,
        active_superpops     : activeRecords.collect { Map record -> record.superpop }.unique().sort(),
        datasets             : activeRecords.collect { Map record -> [id: record.id, assay: record.assay, superpop: record.superpop, sample_count: record.sample_count] }
    ]

    def summaryOut = new File(task.workDir.toString(), 'manifest_summary.json')
    summaryOut.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(summary)) + '\n'
}