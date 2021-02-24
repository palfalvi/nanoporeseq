process unagi {

  label "long_job"

  container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*.gtf'

  input:
    path genome
    path fastq

  output:
    path "*.gtf", emit: gtf

  script:
    def stranded  = params.ont_stranded ? "--stranded"  : ""
    """
    unagi \
    --input $fastq \
    --genome $genome \
    --output unagi \
    $stranded

    $workDir/scripts/bed2gff.py -q unagi/Splicing_Isoforms.bed -o unagi.gff3
    """
}
