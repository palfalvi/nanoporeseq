process braker2 {

  label "small_job"

  conda "$baseDir/conda-envs/braker2-env.yaml"

  publishDir "${params.outdir}/braker2/", mode: 'copy'

  input:
    path genome
    path junctions
    path gtfs
    path blast_db

  output:
    path "mikado_prepared.gtf", emit: gtf
    path "mikado_prepared.fasta", emit: fasta
    path "prepare.log", emit: log
    path "configuration.yaml", emit: config

  script:

    """
    braker.pl \
    --genome=$genome \
    --bam=$bam \
    --prot_seq=$proteins \
    --prg=gth \
    --gth2traingenes \
    --species = $params.species \
    --softmasking \
    --cores $task.cpus \
    --gff3
    """
}
