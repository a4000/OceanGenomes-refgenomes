process TIARA_TIARA {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tiara:1.0.3' :
        'biocontainers/tiara:1.0.3' }"

    input:
    tuple val(meta), path(fasta)
    val(haplotype)

    output:
    tuple val(meta), path("${prefix}*.{txt,txt.gz}"), emit: classifications
    tuple val(meta), path("log_*.{txt,txt.gz}")     , emit: log
    tuple val(meta), path("*.{fasta,fasta.gz}")     , emit: fasta, optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    tiara -i ${fasta} \
        -o ${prefix}_${haplotype}.txt \
        --threads ${task.cpus} \
        ${args}

    find . -name "*_${fasta}*" -exec sh -c 'file=`basename {}`; mv "\$file" "\${file%%_*}_${prefix}_${haplotype}.fasta"' \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiara: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.out.txt
    touch log_${prefix}.out.txt
    touch bacteria_${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiara: ${VERSION}
    END_VERSIONS
    """
}
