process edta {

  label "small_job"

  conda "$baseDir/conda-envs/bedtools-env.yaml"

  publishDir "${params.outdir}/soft_masking", mode: 'copy'

  input:
    path genome
    path tes

  output:
    path "*_softmasked.fasta", emit: masked

  script:
    """
    bedtools maskfasta -soft -fi $genome -bed $tes -fo ${genome.simpleName}_softmasked.fasta
    """
}
