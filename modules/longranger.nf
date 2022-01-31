process longranger {

  label 'long_job'

  container "docker://biocontainers/longranger:v2.2.2_cv2"


  input:
  tuple val(sample_id), file(reads)
  path reference

  output:
    path "${sample_id}/outs/*.bam", emit: bam
    path "${sample_id}/outs/*.bam.bai", emit: baidx

  script:
    """
    longranger align --id=${sample_id} \
                   --reference=${reference} \
                   --fastqs=./
    """
}
