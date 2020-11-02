process wtdbg {

  cpus 32
  memory '1.5T'
  queue 'MPP,smpl'

  conda "$baseDir/conda-envs/wtdbg-env.yaml"

  publishDir "${params.outdir}/wtdbg", mode: 'copy'

  input:
    path fastq
    val genome_size

  output:
    path "wtdbg2_assembly.ctg.fasta", emit: assembly

  script:
    """
    wtdbg2 \
    -x ont \
    -g ${genome_size} \
    -t ${task.cpus} \
    -i $fastq \
    -fo wtdbg2_assembly

    wtpoa-cns -t ${task.cpus} \
    -i wtdbg2_assembly.ctg.lay.gz \
    -fo wtdbg2_assembly.ctg.fasta
    """
}
