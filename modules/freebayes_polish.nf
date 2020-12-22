process nextpolish {

  label 'assembly'

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
