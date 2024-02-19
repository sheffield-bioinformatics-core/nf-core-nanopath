process READ_CLUSTERING {
    tag "$meta.id"
    label (params.throughput == 'high' ? 'high_sensitivity': params.throughput == 'low' ? 'low_resource' : 'standard')

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43ccf4_4' :
        'biocontainers/grep:3.4--hf43ccf4_4' }"

    input:
    tuple val(meta), path(freqs)

    output:
    tuple val(meta), path("*_hdbscan_output.tsv"),              emit: clusters
    path("*_hdbscan_output.png"),                               emit: plot
    path "versions.yml",                                        emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    umap_hdbscan.py \\
        --infile ${freqs} \\
        --outfile ${prefix}_hdbscan_output \\
        $args
    echo "completed clustering"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d ' ' -f2)
        umap: \$(python -c "import umap; print(umap.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}