process transdecoder {
  label "long_job"

  conda "$baseDir/conda-envs/mikado-env.yaml"
  // container "peegee/nanoporeseq:latest"

  // publishDir "${params.outdir}/mikado/", mode: 'copy', pattern: 'mikado*'

  input:
    path fasta

  output:
    path "*transdecoder.pep", emit: pep
    path "*transdecoder.gff3", emit: gff3
    path "*transdecoder.cds", emit: cds
    path "*transdecoder.bed", emit: bed

  script:

    """
    TransDecoder.LongOrfs -t $fasta
    TransDecoder.Predict -t $fasta
    """
}
