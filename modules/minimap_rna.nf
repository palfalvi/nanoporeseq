process minimap2 {

  label 'long_job'

  conda "$baseDir/conda-envs/minimap-env.yaml"

  // publishDir "${params.outdir}/bwa", mode: 'copy'

  input:
    path fastq
    path assembly

  output:
    path "*.bam", emit: bam
    path "*.bam.bai", emit: baidx

  script:
    """
    minimap2 -ax splice -uf -k14 -t ${task.cpus} $assembly $fastq | samtools sort -@ $task.cpus -O BAM - > long_rna.bam
    samtools index long_rna.bam
    """
}
