process scaffX {
  label "small_job"
  publishDir "${params.outdir}/scaff10x", mode: 'copy'

  #conda "$baseDir/conda-envs/scaff10x-env.yaml"

  input:
    path assembly
    tuple val(sample_id), file(reads)

  output:
    path '*.scaff10x.fa', emit: assembly

  script:
    """
    touch reads.dat && \
    echo q1=${reads[0]} > reads.dat && \
    echo q2=${reads[1]} >> reads.dat

    scaff10x \
        -nodes ${task.cpus} \
        -longread 1 \
        -gap 100 \
        -matrix 2000 \
        -reads 10 \
        -score 20 \
        -link 8 \
       -block 50000 \
       -data reads.dat \
        $assembly \
        ${assembly.simpleName}.scaff10x.fa
    """
}
