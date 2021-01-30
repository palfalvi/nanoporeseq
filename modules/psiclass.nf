process psiclass {

  label "small_job"

  //conda "$baseDir/conda-envs/annotation-env.yaml"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*psiclass.gtf'

  input:
    path genome
    tuple file(bam), file(baidx)

  output:
    path "*psiclass.gtf", emit: gtf

  script:

    def strand = params.orientation ? "--fr" : ""

    """
    psiclass \
    -p $task.cpus \
    -o ${bam.simpleName}_psiclass \
    -b $bam

    mv ${bam.simpleName}_psiclass_vote.gtf ${bam.simpleName}_psiclass.gtf
    """
}
