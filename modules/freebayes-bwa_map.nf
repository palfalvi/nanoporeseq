process freebayes_bwa {

  label 'assembly'

  conda "$baseDir/conda-envs/freebayes-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    path assembly
    tuple val(sample_id), file(reads)

  output:
    path "*.bam", emit: bam
    path "*.bam.bai", emit: baidx
    val avg_depth, emit: avg_depth

  script:
    """
    # bwa index
    bwa index ${assembly}

    # bwa map
    bwa mem -t ${task.cpus} ${assembly} ${reads[0]} ${reads[1]} > ${assembly}.sam

    samtools sort -@ ${task.cpus} -O bam -o ${assembly}.bam -T ${assembly}.tmp ${assembly}.sam && rm ${assembly}.sam

    samtools index ${assembly}.bam

    avg_depth=`samtools depth ${assembly}.bam  |  awk '{sum+=\$3} END { print sum/NR}'`
    """
}
