process arima_qc {

  label "small_job"

  conda "$baseDir/conda-envs/samtools-env.yaml"

  input:
    path bam1
    path bam2
    path genome

  output:
    path "filtered*.bam", emit: bam

  script:

    """
    samtools faidx $genome

    perl $baseDir/modules/two_read_bam_combiner.pl $bam1 $bam2 samtools $params.mapq_filter | samtools view -bS -t ${genome}.faidx - | samtools sort -@ $task.cpus -o ${bam1.simpleName}.merged.bam -

    """
}
