process racon {

  cpus 32
  memory '1.5T'
  queue 'MPP,smpl'

  conda "$baseDir/conda-envs/racon-env.yaml"

  // publishDir "${params.outdir}/racon", mode: 'copy'

  input:
    path fastq
    path overlap
    path assembly

  output:
    path "*_racon.fasta", emit: assembly

  script:
    """
    racon \
    --threads ${task.cpus} \
    -u \
    $fastq \
    $overlap \
    $assembly > \
    ${fastq.simpleName}_racon.fasta
    """
}
