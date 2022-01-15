process scaffX {
  tag "$assembly"
  label "small_plus"
  publishDir "${params.outdir}/scaff10x", mode: 'copy'

  input:
    path assembly
    tuple val(sample_id), file(reads)

  output:
    path '*.scaff10x.fa', emit: assembly

  script:
    """
    SCAFF10X=`which scaff10x`

    touch reads.dat && \
    echo q1=${reads[0]} > reads.dat && \
    echo q2=${reads[1]} >> reads.dat

    \$SCAFF10X \
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
