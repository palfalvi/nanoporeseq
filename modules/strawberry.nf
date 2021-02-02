process strawberry {

  label "small_job"

  conda "$baseDir/conda-envs/strawberry-env.yaml"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*.gtf'

  input:
    path genome
    tuple file(bam), file(baidx)

  output:
    path "*.gtf", emit: gtf

  script:

    def strand = params.orientation ? "--fr" : ""

    """
    chmod +x $projectDir/scripts/strawberry

    $projectDir/scripts/strawberry \
    --no-quant \
    $strand \
    --num-threads $task.cpus \
    --output-gtf ${bam.simpleName}_strawberry.gtf \
    $bam

    """
}
