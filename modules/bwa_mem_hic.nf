process bwa_mem_hic {

  label 'assembly'

  conda "$baseDir/conda-envs/bwa-samtools-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    tuple val(sample_id), file(reads)
    val num
    path assembly
    path index

  output:
    path "*.bam", emit: bam
    path "*.bam.bai", emit: baidx

  script:

    def read = num == 0 ? "${reads[0]}" : "${reads[1]}"

    """
    bwa mem -t ${task.cpus} ${assembly} $read | samtools sort -@ ${task.cpus} -O BAM - > ${assembly.simpleName}.bam

    samtools index ${assembly.simpleName}.bam
    """
}
