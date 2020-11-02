process masurca {

  cpus 32
  memory '1.5T'
  queue 'MPP,smpl'

  conda "$baseDir/conda-envs/masurca-env.yaml"

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
