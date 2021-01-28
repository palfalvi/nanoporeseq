process star_idx {
  tag "$genome"
  label 'small_plus'
  cpus "$params.cpus"

  conda "$baseDir/conda-envs/star-env.yaml"

  input:
    path genome
  output:
    path "${genome.simpleName}_idx"
  script:
    """
		mkdir ${genome.simpleName}_idx
    STAR --runMode genomeGenerate \
    --runThreadN ${task.cpus} \
    --genomeDir ${genome.simpleName}_idx \
    --genomeFastaFiles $genome
    """
}
