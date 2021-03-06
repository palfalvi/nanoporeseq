process strawberry {

  label "small_job"

  conda "$baseDir/conda-envs/strawberry-env.yaml"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*.gtf'

  input:
    path genome
    path bam

  output:
    path "*.gtf", emit: gtf

  script:

    def strand = params.orientation ? "--fr" : ""

    """
    chmod +x $projectDir/scripts/strawberry

    $projectDir/scripts/strawberry \
    $strand \
    --num-threads $task.cpus \
    --output-gtf ${bam.simpleName}_strawberry.gtf \
    $bam

    """
}
