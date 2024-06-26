process BUSCO_GENERATEPLOT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'ezlabgva/busco:v5.7.1_cv1'

    input:
    tuple val(meta), path(short_summary_txt, stageAs: 'busco/*')
    val(stage)

    output:
    tuple val(meta), path('*.png'), emit: png
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: 'busco_figure'
    """
    generate_plot.py \\
        $args \\
        -wd busco

    mv ./busco/busco_figure.png ${meta.id}_${prefix}_${stage}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix   ?: 'busco_figure'
    """
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
