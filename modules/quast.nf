process quast {

  cpus 10
  memory '40G'
  queue 'CDE'

  conda "$baseDir/conda-envs/genome-qc.yaml"

  publishDir "${params.outdir}/quast", mode: 'copy'

  input:
    path genomes // Can be a list of genomes

  output:
    path "quast_out/report.tsv", emit: summary

  script:
    """
    quast.py \
    --large \
    --eukaryote \
    --k-mer-stats \
    --threads ${task.cpus} \
    -o quast_out \
    --no-plots \
    --no-icarus \
    --no-html \
    $genomes
    """
}
