process arima_filter {

  label "small_job"

  conda "$baseDir/conda-envs/samtools-env.yaml"

  input:
    path bam

  output:
    path "*filtered.bam", emit: bam

  script:

    """
    samtools view -h $bam | perl $baseDir/scripts/filter_five_end.pl | samtools view -Sb - | samtools sort -@ $task.cpus -o ${bam1.simpleName}.filtered.bam -
    """
}
