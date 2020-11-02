process nextdenovo {

  cpus 32
  memory '1.5T'
  queue 'MPP,smpl'

  //conda "$baseDir/conda-envs/nextdenovo-env.yaml"

  publishDir "${params.outdir}/nextdenovo", mode: 'copy'

  input:
    path nextdenovo
    path cfg

  output:
    path "01_rundir/03.ctg_graph/nd.asm.fasta", emit: assembly
    path "*.gfa", emit: gfa

  script:
    """
    $nextdenovo $cfg
    """
}
