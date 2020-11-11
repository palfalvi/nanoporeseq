process nanolyse {
  label "long_job"

  conda "$baseDir/conda-envs/ont-cleanup.yaml"

  publishDir "${params.outdir}/filtered_reads", mode: 'copy'

  input:
    path reads // A list of fastq files

  output:
    path "*filtered.fastq.gz", emit: filtered
    path "*filter.log", emit: filter_log
    path "*lyse.log", emit: lyse_log

  script:
    def qual      = params.cu_qual     ? "--quality ${params.cu_qual}"       : ""
    def length    = params.cu_length   ? "--length ${params.cu_length}"     : ""
    def headcrop  = params.cu_headcrop ? "--headcrop ${params.cu_headcrop}" : ""

    """
    gunzip -c $reads | NanoLyse --reference $projectDir/conf/ont_control_dna.fasta --logfile ${reads.simpleName}_lyse.log | NanoFilt $qual $length $headcrop --logfile ${reads.simpleName}_filter.log | gzip > ${reads.simpleName}_filtered.fastq.gz
    """
}
