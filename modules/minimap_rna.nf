process minimap_rna {
  tag "$sample_id"
  label 'long_job'

  conda "$baseDir/conda-envs/minimap-env.yaml"

  // publishDir "${params.outdir}/bwa", mode: 'copy'

  input:
    tuple val(sample_id), file(reads)
    path assembly

  output:
    path "*.bam", emit: bam

  script:
    """
    minimap2 -ax splice -uf -k14 -t ${task.cpus} $reads $fastq | samtools sort -@ $task.cpus -O BAM - > ${sample_id}_minimap2.bam
    """
}
