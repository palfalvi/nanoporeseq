process flye {

  label 'assembly'

  conda "$baseDir/conda-envs/flye-env.yaml"

  publishDir "${params.outdir}/flye", mode: 'copy'

  input:
    path fastq
    val genome_size

  output:
    path "flye_out/assembly.fasta", emit: assembly
    path "flye_out/assembly_graph.gfa", emit: gfa

  script:
    """
    flye \
    --nano-raw ${fastq} \
    --threads ${task.cpus} \
    --genome_size ${genome_size} \
    --asm_coverage 50 \
    --out_dir flye_out
    """
}
