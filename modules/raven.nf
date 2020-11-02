process raven {

  cpus 32
  memory '1.5T'
  queue 'MPP'

  conda "$baseDir/conda-envs/raven-env.yaml"

  publishDir "${params.outdir}/raven", mode: 'copy'

  input:
    path fastq

  output:
    path "*.fasta", emit: assembly
    path "*.gfa", emit: gfa

  script:
    """
    raven \
    --threads $task.cpus \
    --graphical-fragment-assembly raven_assembly.gfa \
    $fastq
    """
}
