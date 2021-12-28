process scaff10x {
  label "small_job"

  //conda "$baseDir/conda-envs/scaff10x-env.yaml"

  publishDir "${params.outdir}/scaff10x", mode: 'copy'

  input:
    path scaff10x
    path assembly
    tuple val(sample_id), file(reads)

  output:
    path '*.scaff10x.fa', emit: assembly

  script:
    """
    export PATH=$scaff10x:\$PATH

    echo q1=${reads[1]} >> reads.dat
    echo q2=${reads[2]} >> reads.dat

    $scaff10x/scaff10x \
        -nodes ${task.cpus} \
        -longread 1 \
        -gap 100 \
        -matrix 2000 \
        -reads 10 \
        -score 20 \
        -edge 50000 \
        -link 8 \
        -block 50000 \
        -data reads.dat \
        $assembly \
        ${assembly.simpleName}.scaff10x.fa
    """
}
