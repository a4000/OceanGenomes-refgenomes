process TAR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(dir)
    val(suffix)

    output:
    tuple val(meta), path("*_*"), emit: tar_file
    path  "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tar \\
        $args \\
        ${prefix}_${suffix} \\
        $dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$( tar --version | head -n 1 | sed 's/tar (GNU tar) //g' )
    END_VERSIONS
    """
}
