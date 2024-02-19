process KMER_FREQS {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43ccf4_4' :
        'biocontainers/grep:3.4--hf43ccf4_4' }"

    input:
    tuple val(meta), path(reads)

    output:
    path "*_freqs.txt",                          emit: freqs
    tuple val(meta), path(reads),                emit: reads
    path "versions.yml",                         emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kmer_freq.py -r $reads > ${prefix}_freqs.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version | head -n1 | cut -f 4 -d ' ')
        gzip: \$(gzip --version | head -n1 | cut -f 2 -d ' ')
    END_VERSIONS
    """
}