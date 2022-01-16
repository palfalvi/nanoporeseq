process masurca {

  label 'assembly'

  container 'https://depot.galaxyproject.org/singularity/masurca:4.0.7--pl5262h86ccdc5_0'

  publishDir "${params.outdir}/masurca", mode: 'copy'

  input:
    file masurca_file

  output:
    file "*final.genome.scf.fasta", emit: assembly

  script:
    """
    masurca $masurca_file

    ./assemble.sh

    cp CA.mr*/final.genome.scf.fasta final.genome.scf.fasta
    """
}
