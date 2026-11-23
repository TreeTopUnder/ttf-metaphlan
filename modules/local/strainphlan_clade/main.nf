process STRAINPHLAN_CLADE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::metaphlan=4.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaphlan:4.1.1--pyhdfd78af_0' :
        'biocontainers/metaphlan:4.1.1--pyhdfd78af_0' }"

    input:
    val meta
    path markers, stageAs: 'consensus_markers/*'
    path metaphlan_db

    output:
    path "${meta.id}"  , emit: clade_dir
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--phylophlan_mode fast'
    def clade = meta.id
    """
    mkdir -p db_markers ${clade}

    DB_PKL=\$(find -L "${metaphlan_db}" -name "*.pkl" | head -n 1)
    if [ -z "\$DB_PKL" ]; then
        echo "ERROR: No MetaPhlAn .pkl database found under ${metaphlan_db}" >&2
        exit 1
    fi

    write_status() {
        printf '%s\\nResearch use only. Not for clinical or diagnostic use.\\n' "\$1" > ${clade}/status.txt
    }

    set +e
    extract_markers.py \\
        -c ${clade} \\
        -o db_markers \\
        -d "\$DB_PKL"
    extract_status=\$?
    set -e

    if [ \$extract_status -ne 0 ] || [ ! -s "db_markers/${clade}.fna" ]; then
        write_status "extract_markers.py did not produce db_markers/${clade}.fna (exit \${extract_status}). Skipping StrainPhlAn tree for ${clade}."
    else
        shopt -s nullglob
        marker_files=( consensus_markers/* )
        if [ \${#marker_files[@]} -eq 0 ]; then
            write_status "No consensus marker files staged for ${clade}. Skipping StrainPhlAn tree."
        else
            set +e
            strainphlan \\
                -s consensus_markers/* \\
                -m db_markers/${clade}.fna \\
                -o ${clade} \\
                -c ${clade} \\
                -d "\$DB_PKL" \\
                ${args} \\
                --nprocs ${task.cpus}
            strain_status=\$?
            set -e
            if [ \$strain_status -ne 0 ]; then
                write_status "StrainPhlAn exited \${strain_status} for ${clade} (too few markers/samples after internal QC is the usual cause). Pipeline continues without a tree."
            else
                write_status "ok"
                if ls ${clade}/RAxML_bestTree.*.tre >/dev/null 2>&1; then
                    echo "tree=present" >> ${clade}/status.txt
                fi
            fi
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: \$(metaphlan --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    def clade = meta.id
    """
    mkdir -p ${clade}
    touch ${clade}/RAxML_bestTree.${clade}.StrainPhlAn4.tre
    touch ${clade}/${clade}.StrainPhlAn4_concatenated.aln
    touch ${clade}/${clade}.info
    printf "ok\\nResearch use only. Not for clinical or diagnostic use.\\ntree=present\\n" > ${clade}/status.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        metaphlan: 4.1.1
    END_VERSIONS
    """
}
