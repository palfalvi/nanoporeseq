process bam_merge {

  conda "$baseDir/conda-envs/star-env.yaml"
  //conda "$baseDir/conda-envs/samtools-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    path bams
    val name

  output:
    tuple file("*.merged.bam"), file("*.merged.bam.bai"), emit: bam
    //path "*.merged.bam", emit: bam
    //path "*.merged.bam.bai", emit: baidx

  script:
    """
    samtools merge ${name}.merged.bam $bams

    samtools index ${name}.merged.bam
    """
}
