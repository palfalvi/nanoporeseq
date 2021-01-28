process hisat2_idx {
  tag "$genome"
  label 'small_plus'

  conda "$baseDir/conda-envs/hisat-env.yaml"

  input:
    path genome
  output:
    path "${genome.simpleName}*.ht2*"
  script:
    """
    hisat2-build -p ${task.cpus} $genome ${genome.simpleName}
    """
}
