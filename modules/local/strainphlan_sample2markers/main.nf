process STRAINPHLAN_SAMPLE2MARKERS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::metaphlan=4.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaphlan:4.1.1--pyhdfd78af_0' :
        'biocontainers/metaphlan:4.1.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sam)
    path metaphlan_db

    output:
    tuple val(meta), path("*.json.bz2"), emit: markers, optional: true
    tuple val(meta), path("*.pkl")     , emit: pkl    , optional: true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    DB_PKL=\$(find -L "${metaphlan_db}" -name "*.pkl" | head -n 1)
    if [ -z "\$DB_PKL" ]; then
        echo "ERROR: No MetaPhlAn .pkl database found under ${metaphlan_db}" >&2
        exit 1
    fi

    sample2markers.py \\
        -i ${sam} \\
        -o . \\
        -n ${task.cpus} \\
        -d "\$DB_PKL" \\
        ${args}

    # MetaPhlAn 4.1.1 names consensus markers from the SAM basename.
    if ls *.json.bz2 >/dev/null 2>&1; then
        for marker in *.json.bz2; do
            if [ "\$marker" != "${prefix}.json.bz2" ]; then
                cp -f "\$marker" "${prefix}.json.bz2" || true
            fi
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.json.bz2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: 4.1.1
    END_VERSIONS
    """
}
