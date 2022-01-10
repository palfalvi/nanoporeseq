process debarcodeX {
  label "small_job"

  conda "$baseDir/conda-envs/scaff10x-env.yaml"

  publishDir "${params.outdir}/scaff10x", mode: 'copy'

  input:
    val scaffX
    tuple val(sample_id), file(reads)

  output:
    tuple val(sample_id), file("*BC*.fastq"), emit: fastq

  script:
    """
    git clone  https://github.com/wtsi-hpag/Scaff10X.git
    cd Scaff10X
    ./install.sh
    cd ..

    echo q1=${reads[0]} >> reads.dat
    echo q2=${reads[1]} >> reads.dat

    ./Scaff10X/src/scaff_reads \
        reads.dat \
        ${sample_id}_BC_R1.fastq \
        ${sample_id}_BC_R2.fastq
    """
}