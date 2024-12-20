process CAT_SCAFFOLDS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*_combined_scaffold.fa"), emit: cat_file
    tuple val(meta), path("*H1.scaffold_1.fa")     , emit: hap1_scaffold
    tuple val(meta), path("*H2.scaffold_2.fa")     , emit: hap2_scaffold
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sed 's/scaffold/H1.scaffold/g' *hap1_filtered_scaffolds.fa > ${prefix}_H1.scaffold_1.fa
    sed 's/scaffold/H2.scaffold/g' *hap2_filtered_scaffolds.fa > ${prefix}_H2.scaffold_2.fa

    cat ${prefix}_H1.scaffold_1.fa ${prefix}_H2.scaffold_2.fa > "${prefix}_combined_scaffold.fa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}
