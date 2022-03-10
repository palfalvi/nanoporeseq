process braker2 {

  label "long_job"

  conda "$baseDir/conda-envs/braker2-env.yaml"
  //container "https://depot.galaxyproject.org/singularity/braker2:2.1.6--hdfd78af_5"

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
    def augustus_con = params.augustus_conf_path ? "--AUGUSTUS_CONFIG_PATH=${params.augustus_conf_path}" : ""
    """

    workdir=`pwd`

    braker.pl \
    --genome=$genome \
    $mapping \
    $protein \
    $mark \
    $sp \
    --GENEMARK_PATH=/home/peegee/bin/ \
    --PROTHINT_PATH=/home/peegee/bin/ProtHint-2.6.0/bin/ \
    --CDBTOOLS_PATH=/home/peegee/miniconda3/bin/ \
    --AUGUSTUS_SCRIPTS_PATH=/home/peegee/miniconda3/ \
    --softmasking \
    --cores $task.cpus \
    $augustus_con
    """
}
