process cufflinks {

  label "small_job"

  conda "$baseDir/conda-envs/cufflinks-env.yaml"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*.gtf'

  input:
    path genome
    tuple file(bam), file(baidx)

  output:
    path "*.gtf", emit: gtf

  script:

    """
    cufflinks \
    --num-threads $task.cpus \
    --output-dir /path/to/outputDirectory ${bam.simpleName} \
    $bam

    mv ./${bam.simpleName}/transcripts.gtf ./${bam.simpleName}_cufflinks.gtf
    """
}
