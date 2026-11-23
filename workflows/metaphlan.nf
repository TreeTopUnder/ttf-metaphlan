include { FASTP } from '../modules/nf-core/fastp/main'
include { BOWTIE2_ALIGN } from '../modules/nf-core/bowtie2/align/main'
include { UNTAR as UNTAR_HOST_INDEX } from '../modules/nf-core/untar/main'
include { UNTAR as UNTAR_METAPHLAN_DB } from '../modules/nf-core/untar/main'
include { METAPHLAN_METAPHLAN } from '../modules/nf-core/metaphlan/metaphlan/main'
include { METAPHLAN_MERGEMETAPHLANTABLES } from '../modules/nf-core/metaphlan/mergemetaphlantables/main'
include { KRONA_KTIMPORTTEXT as KRONA_PER_SAMPLE } from '../modules/nf-core/krona/ktimporttext/main'
include { KRONA_KTIMPORTTEXT as KRONA_COMBINED } from '../modules/nf-core/krona/ktimporttext/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { METAPHLAN_TIDY_TABLES } from '../modules/local/metaphlan_tidy_tables/main'
include { METAPHLAN_TO_KRONA } from '../modules/local/metaphlan_to_krona/main'
include { BUILD_PHYLOSEQ } from '../modules/local/build_phyloseq/main'

def isArchivePath(path) {
    // Unpacked dirs (local or s3://, trailing slash optional) skip UNTAR.
    return RefPaths.isArchive(path)
}

workflow TTF_METAPHLAN {
    take:
    ch_input
    ch_metadata
    ch_metaphlan_db_archive
    ch_metaphlan_db_dir
    ch_host_index_archive
    ch_host_index_dir

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_bowtie2_logs  = Channel.empty()

    FASTP(
        ch_input,
        [],
        false,
        false,
        false
    )
    ch_versions      = ch_versions.mix(FASTP.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] })
    ch_reads         = FASTP.out.reads

    if (isArchivePath(params.metaphlan_db)) {
        UNTAR_METAPHLAN_DB(ch_metaphlan_db_archive)
        ch_metaphlan_db = UNTAR_METAPHLAN_DB.out.untar.map { _meta, dir -> dir }
        ch_versions     = ch_versions.mix(UNTAR_METAPHLAN_DB.out.versions)
    } else {
        ch_metaphlan_db = ch_metaphlan_db_dir
    }

    if (!params.skip_host_removal) {
        if (isArchivePath(params.host_index)) {
            UNTAR_HOST_INDEX(ch_host_index_archive)
            ch_host_index = UNTAR_HOST_INDEX.out.untar
            ch_versions   = ch_versions.mix(UNTAR_HOST_INDEX.out.versions)
        } else {
            ch_host_index = ch_host_index_dir
        }

        BOWTIE2_ALIGN(
            ch_reads,
            ch_host_index,
            [[:], []],
            true,
            false
        )
        ch_versions      = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect { it[1] })
        ch_bowtie2_logs  = BOWTIE2_ALIGN.out.log.map { _meta, log -> log }
        ch_reads         = BOWTIE2_ALIGN.out.fastq
    }

    ch_metaphlan_db_value = ch_metaphlan_db.collect().map { dirs -> dirs[0] }

    METAPHLAN_METAPHLAN(
        ch_reads,
        ch_metaphlan_db_value,
        false
    )
    ch_versions      = ch_versions.mix(METAPHLAN_METAPHLAN.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(METAPHLAN_METAPHLAN.out.profile.collect { it[1] })

    ch_profiles_merged = METAPHLAN_METAPHLAN.out.profile
        .map { _meta, profile -> profile }
        .collect()
        .map { profiles -> [[id: 'merged'], profiles] }

    METAPHLAN_MERGEMETAPHLANTABLES(ch_profiles_merged)
    ch_versions = ch_versions.mix(METAPHLAN_MERGEMETAPHLANTABLES.out.versions)

    ch_logs_for_tidy = ch_bowtie2_logs.collect().ifEmpty([file("${projectDir}/assets/NO_BT2_LOG")])
    METAPHLAN_TIDY_TABLES(
        METAPHLAN_MERGEMETAPHLANTABLES.out.txt,
        ch_logs_for_tidy
    )
    ch_versions = ch_versions.mix(METAPHLAN_TIDY_TABLES.out.versions)

    if (!params.skip_phyloseq) {
        BUILD_PHYLOSEQ(
            METAPHLAN_TIDY_TABLES.out.species,
            METAPHLAN_TIDY_TABLES.out.sgb,
            ch_metadata,
            ch_metaphlan_db_value
        )
        ch_versions = ch_versions.mix(BUILD_PHYLOSEQ.out.versions)
    }

    if (!params.skip_krona) {
        METAPHLAN_TO_KRONA(METAPHLAN_METAPHLAN.out.profile)
        ch_versions = ch_versions.mix(METAPHLAN_TO_KRONA.out.versions)

        KRONA_PER_SAMPLE(METAPHLAN_TO_KRONA.out.txt)
        ch_versions = ch_versions.mix(KRONA_PER_SAMPLE.out.versions)

        ch_krona_combined = METAPHLAN_TO_KRONA.out.txt
            .map { _meta, txt -> txt }
            .collect()
            .map { files -> [[id: 'krona_combined'], files] }
        KRONA_COMBINED(ch_krona_combined)
        ch_versions = ch_versions.mix(KRONA_COMBINED.out.versions)
    }

    ch_collated_versions = ch_versions
        .map { path -> path.text }
        .unique()
        .mix(
            Channel.of(
                """
                Workflow:
                    ttf-metaphlan: ${workflow.manifest.version}
                    nextflow: ${workflow.nextflow.version}
                """.stripIndent()
            )
        )
        .collectFile(
            name: 'software_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            sort: true,
            newLine: true
        )
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    MULTIQC(
        ch_multiqc_files.collect(),
        file("${projectDir}/assets/multiqc_config.yml"),
        [],
        [],
        [],
        []
    )

    emit:
    versions = ch_versions
}
