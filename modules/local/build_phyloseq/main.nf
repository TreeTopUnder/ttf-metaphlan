process BUILD_PHYLOSEQ {
    tag "phyloseq"
    label 'process_low'

    conda "bioconda::bioconductor-phyloseq=1.46.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-phyloseq:1.46.0--r43hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-phyloseq:1.46.0--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(species_table)
    tuple val(meta2), path(sgb_table)
    path metadata
    path metaphlan_db

    output:
    path "metaphlan_phyloseq.rds", emit: rds
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def metadata_arg = metadata ? "--metadata ${metadata}" : ""
    def db_arg = metaphlan_db ? "--db ${metaphlan_db}" : ""
    """
    build_phyloseq.R \\
        --species ${species_table} \\
        --sgb ${sgb_table} \\
        ${metadata_arg} \\
        ${db_arg} \\
        --output metaphlan_phyloseq.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | head -n 1 | sed 's/R version //; s/ .*//')
        phyloseq: \$(Rscript -e 'cat(as.character(packageVersion("phyloseq")))')
    END_VERSIONS
    """

    stub:
    """
    touch metaphlan_phyloseq.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: 4.3.0
        phyloseq: 1.46.0
    END_VERSIONS
    """
}
