process masurca {

  label 'assembly'

  container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/masurca", mode: 'copy'

  input:
    file masurca_file

  output:
    file "final.genome.scf.fasta", emit: assembly

  script:
    """
    masurca $masurca_file

    ./assemble.sh

    cp CA.mr*/final.genome.scf.fasta final.genome.scf.fasta
    """
}
