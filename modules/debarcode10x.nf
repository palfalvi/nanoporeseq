process debarcodeX {
  label "small_job"

  //conda "$baseDir/conda-envs/scaff10x-env.yaml"

  publishDir "${params.outdir}/scaff10x", mode: 'copy'

  input:
    path scaff10x
    tuple val(sample_id), file(reads)

  output:
    tuple val(sample_id), file("*BC*.fastq"), emit: fastq

  script:
    """
    export $scaff10x:\$PATH

    echo q1=${reads[1]} >> reads.dat
    echo q2=${reads[2]} >> reads.dat

    $scaff10x/scaff_reads \
        reads.dat \
        ${sample_id}_BC_R1.fastq \
        ${sample_id}_BC_R2.fastq
    """
}
