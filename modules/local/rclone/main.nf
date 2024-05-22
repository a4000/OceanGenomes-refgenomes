process RCLONE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "hub.docker.com/rclone/rclone:1.65.2"

    input:
    tuple val(meta), path(files)
    val(dest_path)

    output:
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    files=$files
    if [[ \${#files[@]} -eq 1 ]]; then
        rclone \\
            move \\
            $args \\
            \$files \\
            $dest_path
    else
        for file in "\${files[@]}"; do
            rclone \\
            move \\
            $args \\
            \$file \\
            $dest_path
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rclone: \$( rclone --version | head -n 1 | sed 's/rclone //g' )
    END_VERSIONS
    """
}
