process unagi {

  label "long_job"

  container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*.gtf'

  input:
    path genome
    tuple val(sample_id), file(reads)

  output:
    path "*.gff3", emit: gtf

  script:
    def stranded  = params.ont_stranded ? "--stranded"  : ""
    """
    gunzip -c $reads > ${reads.simpleName}.fastq

    unagi \
    --input ${reads.simpleName}.fastq \
    --genome $genome \
    --output unagi \
    $stranded

    $baseDir/scripts/bed2gff.py -q unagi/Splicing_Isoforms.bed -o unagi_${sample_id}.gff3
    """
}
