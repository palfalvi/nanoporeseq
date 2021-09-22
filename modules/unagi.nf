process unagi {

  label "long_job"

  conda "$baseDir/conda-envs/braker2-env.yaml"
  // container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*.gtf'

  input:
    path genome
    tuple val(sample_id), file(reads)

  output:
    path "unagi/Splicing_Isoforms.bed", emit: gtf

  script:
    def stranded  = params.ont_stranded ? "--stranded"  : ""
    """

    wget https://github.com/iMetOsaka/UNAGI/archive/refs/heads/master.zip
    unzip master.zip
    chmod a+x UNAGI-master/unagi

    if [[ $reads == *.gz ]]
    then
      gunzip -c $reads > ${reads.simpleName}.fastq
    else
      cat $reads > ${reads.simpleName}.fastq
    fi

    UNAGI-master/unagi \
    --input ${reads.simpleName}.fastq \
    --genome $genome \
    --output unagi \
    $stranded

    """
}
