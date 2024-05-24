process RCLONE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/toolhippie/rclone:20240212"

    input:
    tuple val(meta), path(files), val(dest_path)

    output:
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir files_to_move
    cp $files files_to_move
    cd files_to_move

    rclone move \\
        \$(pwd)/. \\
        $dest_path \\
        $args

    cd ../

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rclone: \$( rclone --version | head -n 1 | sed 's/rclone //g' )
    END_VERSIONS
    """
}
