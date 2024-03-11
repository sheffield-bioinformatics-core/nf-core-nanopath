process FULL_CLASSIFICATION {
    tag "$meta.id"+ "_" + "$cluster"

    conda "bioconda::kraken2 rdptools blast"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/grep:3.4--hf43ccf4_4' :
        'docker.io/mbdabrowska1/full-classification:1.0' }"

    input:
    tuple val(meta), path(consensus), path(cluster_log), val(cluster)

    output:
    tuple val(meta), path('*_consensus_classification.csv'),      emit: classification
    tuple val(meta), path('*_classification.log'),                emit: log
    tuple val(meta), path('*_classification_out.tsv'),                 emit: tsv
    path "versions.yml",                                                                       emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cluster_id = "${cluster}"
    def kraken2_db = params.kraken2_db
    def seqmatch_db = params.seqmatch_db
    def seqmatch_accession = params.seqmatch_accession
    def blast_db = params.blast_db

    """
    echo "chosen classification: full"
    echo "classifying with kraken2"
    kraken2 --db ${kraken2_db} --report ${prefix}_${cluster_id}_kraken2_consensus_classification.csv --output ${prefix}_${cluster_id}_kraken2_classification_out.tsv $consensus
    KR_OUT=\$(sed 's/\t/;/g' ${prefix}_${cluster_id}_kraken2_consensus_classification.csv | tr -s ' ' | sed 's/; /;/g' | cut -d ';' -f3,4,5,6 | grep -v '^0' | awk 'BEGIN {FS=";"; OFS=";"} {print \$4, \$3, \$2}')
    
    echo "classifying with seqmatch"
    SequenceMatch seqmatch -k 5 ${seqmatch_db} $consensus | cut -f2,4 | sort | join -t \$'\t' -1 1 -2 1 -o 2.3,2.5,1.2 - ${seqmatch_accession} | sort -k3 -n -r -t '\t' | sed 's/\t/;/g' > ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv
    if [ -s ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv ]; then
        echo "success"
    else
        echo "unclassified;0;0" >> ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv
    fi
    SEQ_OUT=\$(head -n1 ${prefix}_${cluster_id}_seqmatch_consensus_classification.csv)

    echo "classifying with blastn"
    export BLASTDB=\$(dirname ${blast_db})
    blastn -query $consensus -db \$(basename ${blast_db}) -task megablast -dust no -outfmt "10 sscinames staxids evalue length pident bitscore" -evalue 11 -max_hsps 50 -max_target_seqs 5 | sed 's/,/;/g' > ${prefix}_${cluster_id}_blastn_consensus_classification.csv
    if [ -s ${prefix}_${cluster_id}_blastn_consensus_classification.csv ]; then
        echo "success"
    else
        echo "unclassified;0;0" >> ${prefix}_${cluster_id}_blastn_consensus_classification.csv
    fi
    BLAST_OUT=\$(cut -d";" -f1,2,5 ${prefix}_${cluster_id}_blastn_consensus_classification.csv | head -n1)

    cat $cluster_log > ${prefix}_${cluster_id}_classification.log
    echo -n ";" >> ${prefix}_${cluster_id}_classification.log
    echo \$KR_OUT >> ${prefix}_${cluster_id}_classification.log
    echo \$SEQ_OUT >> ${prefix}_${cluster_id}_classification.log
    echo \$BLAST_OUT >> ${prefix}_${cluster_id}_classification.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version | head -n1 | cut -d' ' -f3)
        blastn: \$(blastn -version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """


}