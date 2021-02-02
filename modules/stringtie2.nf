process stringtie2 {

  label "small_job"

  conda "$baseDir/conda-envs/stringtie2-env.yaml"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*.gtf'

  input:
    path genome
    tuple file(bam), file(baidx)
    val extra

  output:
    path "*.gtf", emit: gtf

  script:

    """
    stringtie \
    -p $task.cpus \
    -o ${bam.simpleName}_stringtie2.gtf \
    $extra \
    $bam

    """
}
