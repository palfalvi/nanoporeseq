process braker2 {

  label "long_job"

  conda "$baseDir/conda-envs/braker2-env.yaml"
  //container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: 'braker/braker.gtf'

  input:
    path genome
    path bam
    val mark

  output:
    path "braker/*augustus.hints.gtf", emit: gtf

  script:
    def protein   = params.protein     ? "--prot_seq=${params.protein}"               : ""
    def mapping   = bam!=[]          ? "--bam ${bam}"                               : ""
    def sp        = params.species     ? "--species ${params.species} --useexisting"  : ""
    def augustus_con = "--AUGUSTUS_CONFIG_PATH=${params.augustus_conf_path}"
    """

    workdir=`pwd`

    braker.pl \
    --genome=$genome \
    $mapping \
    $protein \
    $mark \
    $sp \
    --softmasking \
    --cores $task.cpus \
    $augustus_con
    """
}
