process guppy_basecalling {

  label 'long_job'

  publishDir "${params.outdir}/guppy", mode: 'copy'

  input:
    tuple val(sample_id), file(reads)
    path guppy

  output:
    path "$sample_id/guppy"

  script:
    """
    mkdir $sample_id
    mkdir $sample_id/raw

    tar -xvzf $reads --directory $sample_id/raw

    $guppy \
    -i $sample_id/raw \
    --recursive \
    -s $sample_id/guppy \
    --config  $params.config_file \
    --cpu_threads_per_caller 5 \
    --num_callers 2 \
    --records_per_fastq 100000000 \
    --bam_out ${sample_id}_methylation.bam \
    --qscore_filtering \
    --min_qscore 7 \
    --compress_fastq
    """
}
