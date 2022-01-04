process scaffX {
  label "small_job"

  conda "$baseDir/conda-envs/scaff10x-env.yaml"

  publishDir "${params.outdir}/scaff10x", mode: 'copy'

  input:
    path scaffX
    path assembly
    tuple val(sample_id), file(reads)

  output:
    path '*.scaff10x.fa', emit: assembly

  script:
    """
    git clone  https://github.com/wtsi-hpag/Scaff10X.git
    cd Scaff10X
    ./install.sh
    cd ..

    echo q1=${reads[0]} > reads.dat
    echo q2=${reads[1]} >> reads.dat

    ./Scaff10X/src/scaff10x \
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
