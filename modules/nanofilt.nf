process nanolyse {
  label "small_job"

  conda "$baseDir/conda-envs/ont-cleanup.yaml"

  publishDir "${params.outdir}/filtered_reads", mode: 'copy'

  input:
    path reads // A list of fastq files

  output:
    path "*.fastq"

  script:
    def ref       = params.quast_reference  ? "-r ${params.quast_reference}"        : ""
    def features  = params.quast_features   ? "--features ${params.quast_features}" : ""

    """
    nanolyse $reads
    """
}
