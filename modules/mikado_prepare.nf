process mikado_prepare {

  label "small_job"

  conda "$baseDir/conda-envs/annotation-env.yaml"

  publishDir "${params.outdir}/mikado/conf/", mode: 'copy'

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
    mikado configure \
    --list list.txt \
    --reference $genome \
    --mode permissive \
    --scoring plants.yaml  \
    --junctions $junctions \
    -bt $blast_db uniprot_sprot_plants.fasta \
    configuration.yaml

    mikado prepare \
    --json-conf configuration.yaml
    """
}
