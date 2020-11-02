process nextpolish {

  cpus 32
  memory '1.5T'
  queue 'MPP,smpl'

  conda "$baseDir/conda-envs/freebayes-env.yaml"

  // publishDir "${params.outdir}/nextpolish", mode: 'copy'

  input:
    path assembly
    path reads

  output:
    path "*nextpolish.fasta", emit: assembly

  script:
    """

    """
}
