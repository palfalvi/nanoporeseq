process braker2 {

  label "long_job"

  //conda "$baseDir/conda-envs/braker2-env.yaml"
  container "palfalvi/nanoporeseq"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy'

  input:
    path genome
    path bam
    val mark

  output:
    path "*braker.gtf", emit: gtf

  script:
    def protein   = params.protein     ? "--prot_seq=${params.protein}"               : ""
    def mapping   = bam!=null          ? "--bam ${bam}"                               : ""
    def sp        = params.species     ? "--species ${params.species} --useexisting"  : ""
    """
    braker.pl \
    --genome=$genome \
    $mapping \
    $protein \
    --prg=gth \
    $mark \
    --gth2traingenes \
    $sp \
    --softmasking \
    --cores $task.cpus
    """
}
