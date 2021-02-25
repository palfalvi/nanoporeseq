process psiclass {

  label "small_job"

  //conda "$baseDir/conda-envs/psiclass-env.yaml"
  container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*psiclass.gtf'

  input:
    path genome
    tuple file(bam), file(baidx)

  output:
    path "psiclass_vote.gtf", emit: gtf

  script:

    """
    ls -1 *.bam > bamlist.txt


    psiclass \
    -p $task.cpus \
    --lb bamlist.txt
    """
}
