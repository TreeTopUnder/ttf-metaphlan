process METAPHLAN_TIDY_TABLES {
    tag "merged"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(merged_table)
    path bowtie2_logs

    output:
    tuple val(meta), path("merged_species_table.tsv"), emit: species
    tuple val(meta), path("merged_sgb_table.tsv")    , emit: sgb
    path "host_removal_summary.tsv"                  , emit: host_summary, optional: true
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    shopt -s nullglob
    LOG_ARGS=""
    LOGS=( *.bowtie2.log bt2logs/*.bowtie2.log bt2logs/*.log )
    if [ \${#LOGS[@]} -gt 0 ]; then
        LOG_ARGS="--bowtie2-logs \${LOGS[*]} --host-summary-out host_removal_summary.tsv"
    fi

    tidy_metaphlan_tables.py \\
        --merged ${merged_table} \\
        --species-out merged_species_table.tsv \\
        --sgb-out merged_sgb_table.tsv \\
        \$LOG_ARGS

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf "clade_name\\tsample1\\n" > merged_species_table.tsv
    printf "clade_name\\tsample1\\n" > merged_sgb_table.tsv
    shopt -s nullglob
    logs=( *.bowtie2.log )
    if [ \${#logs[@]} -gt 0 ]; then
        printf "sample\\ttotal_reads\\taligned_reads\\tunaligned_reads\\toverall_alignment_rate\\n" > host_removal_summary.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
