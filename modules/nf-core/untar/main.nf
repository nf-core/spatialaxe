process UNTAR {
    tag "${archive}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${prefix}"), emit: untar
    tuple val("${task.process}"), val('untar'), eval('tar --version 2>&1 | head -1 | sed "s/tar (GNU tar) //; s/ Copyright.*//"'), emit: versions_untar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: (meta.id ? "${meta.id}" : archive.baseName.toString().replaceFirst(/\.tar$/, ""))

    """
    mkdir ${prefix}

    ## Ensures --strip-components only applied when top level of tar contents is a directory
    ## If just files or multiple directories, place all in prefix
    if [[ \$(tar -taf ${archive} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
        tar \\
            -C ${prefix} --strip-components 1 \\
            -xavf \\
            ${args} \\
            ${archive} \\
            ${args2}
    else
        tar \\
            -C ${prefix} \\
            -xavf \\
            ${args} \\
            ${archive} \\
            ${args2}
    fi

    """

    stub:
    prefix = task.ext.prefix ?: (meta.id ? "${meta.id}" : archive.toString().replaceFirst(/\.[^\.]+(.gz)?$/, ""))
    """
    mkdir ${prefix}
    ## Emit a valid empty file. nf-test's snapshot md5 path decompresses .gz files,
    ## so a 0-byte .gz (from `touch`) throws EOF in GZIPInputStream and the snapshot
    ## falls back to dumping non-deterministic File metadata (freeSpace, etc).
    ## `: | gzip -n` produces a 20-byte deterministic empty gzip that md5s consistently.
    emit_stub() {
        local f="\$1"
        if [[ "\$f" == *.gz ]]; then
            : | gzip -n > "\$f"
        else
            touch "\$f"
        fi
    }

    ## Dry-run untaring the archive to get the files and place all in prefix
    if [[ \$(tar -taf ${archive} | grep -o -P "^.*?\\/" | uniq | wc -l) -eq 1 ]]; then
        ## Single top-level dir in the archive: mirror `--strip-components 1` from the real
        ## extraction so the stubbed file layout matches what a real run would produce.
        for i in `tar -tf ${archive}`;
        do
            stripped=\${i#*/}
            [ -z "\$stripped" ] && continue
            if [[ \$(echo "\${i}" | grep -E "/\$") == "" ]];
            then
                mkdir -p ${prefix}/\$(dirname "\$stripped")
                emit_stub ${prefix}/\$stripped
            else
                mkdir -p ${prefix}/\$stripped
            fi
        done
    else
        for i in `tar -tf ${archive}`;
        do
            if [[ \$(echo "\${i}" | grep -E "/\$") == "" ]];
            then
                emit_stub ${prefix}/\${i}
            else
                mkdir -p ${prefix}/\${i}
            fi
        done
    fi
    """
}
