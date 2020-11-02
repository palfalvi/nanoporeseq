process busco {

  cpus 10
  memory '40G'
  queue 'CDE'

  conda "$baseDir/conda-envs/genome-qc.yaml"

  publishDir "${params.outdir}/busco", mode: 'copy'

  input:
    path genome
    val lineage
    val mode

  output:
    path "busco_${genome.simpleName}_${lineage}/short_summary*", emit: summary

  script:
    """
    busco \
    --in $genome \
    --lineage_dataset $lineage \
    --out busco_${genome.simpleName}_${lineage} \
    --mode $mode \
    --long \
    --cpu ${task.cpus}
    """
}
