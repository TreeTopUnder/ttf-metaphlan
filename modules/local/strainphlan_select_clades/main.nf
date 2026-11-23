process STRAINPHLAN_SELECT_CLADES {
    tag "clades"
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(sgb_table)
    val clades
    val top_n
    val min_samples

    output:
    path "strainphlan_summary.tsv", emit: summary
    path "selected_clades.txt"    , emit: clade_list
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def clades_arg = clades ? clades.toString() : ''
    """
    select_strainphlan_clades.py \\
        --sgb-table ${sgb_table} \\
        --clades "${clades_arg}" \\
        --top-n ${top_n} \\
        --min-samples ${min_samples} \\
        --summary-out strainphlan_summary.tsv \\
        --clades-out selected_clades.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def clades_arg = clades ? clades.toString() : ''
    """
    select_strainphlan_clades.py \\
        --sgb-table ${sgb_table} \\
        --clades "${clades_arg}" \\
        --top-n ${top_n} \\
        --min-samples ${min_samples} \\
        --summary-out strainphlan_summary.tsv \\
        --clades-out selected_clades.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11
    END_VERSIONS
    """
}
