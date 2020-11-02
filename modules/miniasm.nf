process miniasm {

  cpus 32
  memory '1.5T'
  queue 'MPP,smpl'

  conda "$baseDir/conda-envs/miniasm-env.yaml"

  publishDir "${params.outdir}/miniasm", mode: 'copy'

  input:
    path fastq

  output:
    path "miniasm_assembly.fasta", emit: assembly
    path "miniasm_assembly.gfa", emit: gfa

  script:
    """
    minimap2/minimap2 -x ava-ont -t ${task.cpus} ${fastq} ${fastq} | gzip -1 > reads.paf.gz

    miniasm/miniasm -f ${fastq} reads.paf.gz > miniasm_assembly.gfa

    awk '$1 ~/S/ {print ">"$2"\n"$3}' miniasm_assembly.gfa > miniasm_assembly.fasta
    """
}