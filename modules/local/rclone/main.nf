process RCLONE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "quay.io/toolhippie/rclone:20240212"

    input:
    tuple val(meta), val(step), path(files), val(dest_path)

    output:
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${step}"
    """
    mkdir files_to_move
    cp -rL $files files_to_move
    cd files_to_move
    rclone size . > ${prefix}_rclone_size_before_transfer.txt

    rclone move \\
        \$(pwd)/. \\
        $dest_path \\
        $args

    cd ../
    rm -r files_to_move

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rclone: \$( rclone --version | head -n 1 | sed 's/rclone //g' )
    END_VERSIONS
    """
}
