process bwa_index {

  label 'large_mem'

  conda "$baseDir/conda-envs/bwa-samtools-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    path assembly

  output:
    path "${assembly}.*", emit: index

  script:
    """
    bwa index ${assembly}
    """
}
