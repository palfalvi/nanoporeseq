process bwa_mem_hic {

  label 'assembly'

  conda "$baseDir/conda-envs/bwa-samtools-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    tuple val(sample_id), file(reads)
    value num
    path assembly
    path index

  output:
    path "*.bam", emit: bam
    path "*.bam.bai", emit: baidx

  script:
    """
    bwa mem -t ${task.cpus} ${assembly} ${reads[$num]} | samtools sort -@ ${task.cpus} -O BAM - > ${assembly.simpleName}.bam

    samtools index ${assembly.simpleName}.bam
    """
}
