#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TTF_METAPHLAN } from './workflows/metaphlan'

def helpMessage() {
    log.info(
        """
        ttf-metaphlan ${workflow.manifest.version}
        Shotgun Taxonomic Profiling (MetaPhlAn 4) for TreeTopFlow

        Research use only. Not for clinical or diagnostic use.

        Required:
          --input             Samplesheet CSV with columns sample,fastq_1,fastq_2
          --outdir            Output directory

        Optional:
          --metadata          Sample metadata TSV (first column ID or sample)
          --metaphlan_db      MetaPhlAn 4 DB dir or .tar/.tar.gz (default: s3://codonyat/refdb/public/metaphlan4/mpa_vJun23_CHOCOPhlAnSGB_202403)
          --host_index        Host Bowtie2 index dir or .tar/.tar.gz (default: s3://codonyat/refdb/public/host/GRCh38_noalt_as)
          --skip_host_removal Skip host read removal (default: false)
          --analysis_type     rel_ab_w_read_stats | rel_ab
          --add_viruses       Include viral clades (default: false)
          --unclassified_estimation  Estimate unclassified fraction (default: true)
          --min_read_len      fastp --length_required (default: 70)
          --skip_krona        Skip Krona plots (default: false)
          --skip_phyloseq     Skip phyloseq RDS (default: false)
          --run_strainphlan   Run StrainPhlAn 4 strain-level profiling (default: false)
          --strainphlan_clades        Optional comma-separated SGB IDs (e.g. t__SGB1877,t__SGB15316)
          --strainphlan_top_n         Max auto-selected SGB clades (default: 5)
          --strainphlan_min_samples   Min samples carrying an SGB for tree building (default: 4)
          --strainphlan_phylophlan_mode  PhyloPhlAn mode: fast | accurate (default: fast)
        """.stripIndent()
    )
}

def isArchive(path) {
    return RefPaths.isArchive(path)
}

def asRefFile(path) {
    def value = RefPaths.normalize(path)
    if (!value) {
        return null
    }
    // Unpacked dirs (local or s3://, trailing slash optional) stage as directories.
    return isArchive(value) ? file(value) : file(value, type: 'dir')
}

def resolveInputPath(raw) {
    if (raw == null) {
        return null
    }
    def value = raw.toString().trim()
    if (!value) {
        return null
    }
    def candidate = file(value)
    if (candidate.exists()) {
        return candidate
    }
    def fromProject = file("${workflow.projectDir}/${value}")
    if (fromProject.exists()) {
        return fromProject
    }
    return candidate
}

def createFastqChannel(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.sample?.toString()?.trim()
    if (!meta.id) {
        error "Samplesheet rows must include a non-empty 'sample' value"
    }

    def fastq_1 = resolveInputPath(row.fastq_1)
    if (!fastq_1) {
        error "Sample ${meta.id} is missing fastq_1"
    }

    def fastq_2 = resolveInputPath(row.containsKey('fastq_2') ? row.fastq_2 : null)
    meta.single_end = fastq_2 == null
    if (meta.single_end) {
        return [meta, [fastq_1]]
    }
    return [meta, [fastq_1, fastq_2]]
}

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    if (!params.input) {
        error "Missing required parameter --input (samplesheet CSV)"
    }
    if (!params.outdir) {
        error "Missing required parameter --outdir"
    }
    if (!params.metaphlan_db) {
        error "Missing required parameter --metaphlan_db"
    }
    if (!params.skip_host_removal && !params.host_index) {
        error "Missing required parameter --host_index (required unless --skip_host_removal)"
    }
    if (!(params.analysis_type in ['rel_ab_w_read_stats', 'rel_ab'])) {
        error "Parameter --analysis_type must be one of: rel_ab_w_read_stats, rel_ab"
    }
    if (!(params.strainphlan_phylophlan_mode in ['fast', 'accurate'])) {
        error "Parameter --strainphlan_phylophlan_mode must be one of: fast, accurate"
    }

    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> createFastqChannel(row) }

    ch_metadata = params.metadata ? Channel.value(file(params.metadata, checkIfExists: true)) : Channel.value([])

    def metaphlan_db_file = asRefFile(params.metaphlan_db)
    ch_metaphlan_db_archive = isArchive(params.metaphlan_db) ? Channel.of([[id: 'metaphlan_db'], metaphlan_db_file]) : Channel.empty()
    ch_metaphlan_db_dir     = isArchive(params.metaphlan_db) ? Channel.empty() : Channel.value(metaphlan_db_file)

    def host_index_file = params.host_index ? asRefFile(params.host_index) : null
    ch_host_index_archive = (!params.skip_host_removal && host_index_file && isArchive(params.host_index)) ? Channel.of([[id: 'host_index'], host_index_file]) : Channel.empty()
    ch_host_index_dir     = (!params.skip_host_removal && host_index_file && !isArchive(params.host_index)) ? Channel.value([[id: 'host_index'], host_index_file]) : Channel.empty()

    TTF_METAPHLAN(
        ch_input,
        ch_metadata,
        ch_metaphlan_db_archive,
        ch_metaphlan_db_dir,
        ch_host_index_archive,
        ch_host_index_dir
    )
}

workflow.onComplete {
    log.info(workflow.success ? "ttf-metaphlan completed successfully (research use only)" : "ttf-metaphlan failed")
}
