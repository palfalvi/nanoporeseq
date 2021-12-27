process medaka {

  label 'assembly'

  conda "$baseDir/conda-envs/medaka-env.yaml"

  publishDir "${params.outdir}/medaka", mode: 'copy'

  input:
    path fastq
    path assembly

  output:
    path "medaka_polish/consensus.fasta", emit: assembly

  script:
    def medaka_model = params.medaka_model ? "-m ${params.medaka_model}" : ""
    """
    medaka_consensus \
    -i $fastq \
    -d $assembly \
    -o medaka_polish \
    -t $task.cpus \
    $medaka_model
    """
}
