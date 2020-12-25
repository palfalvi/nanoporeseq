process pychopper {
  label "long_job"

  tag "$reads"

  conda "$baseDir/conda-envs/pychopper.yaml"

  publishDir "${params.outdir}/filtered_reads", mode: 'copy'

  input:
    path reads // A list of fastq files

  output:
    path "*filtered.fastq", emit: filtered
    path "*unclassified.fastq", emit: unclassified
    path "*rescued.fastq", emit: rescued
    path "*report.pdf", emit: report

  script:

    def qual       = params.cu_qual             ? "-Q ${params.cu_qual}"         : ""
    def length     = params.cu_length           ? "-z ${params.cu_length}"       : ""
    def phmm_file  = params.cu_phmm_file        ? "-g ${params.cu_phmm_file}"    : ""
    def adapter    = params.cu_adapters         ? "-b ${params.cu_adapters}"     : ""
    def config     = params.cu_primer_config    ? "-c ${params.cu_primer_config}": ""
    def method     = params.cu_method           ? "-m ${params.cu_method}"       : ""

    """
    cdna_classifier.py -t $task.cpus -r ${reads.simpleName}_report.pdf $method $adapter $config $phmm_file $qual $length -u ${reads.simpleName}_unclassified.fastq -w ${reads.simpleName}_rescued.fastq $reads ${reads.simpleName}_filtered.fastq
    """
}
