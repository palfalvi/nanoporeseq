process mikado_pick {
  label "long_job"

  conda "$baseDir/conda-envs/mikado-env.yaml"
  // container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/mikado/", mode: 'copy', pattern: 'mikado*'

  input:
    path genome
    path prepare
    path serialise
    path scoring
    // path junction

  output:
    path "mikado.loci.gff3", emit: loci
    path "mikado.subloci.gff3", emit: subloci
    path "mikado.loci.metrics.tsv", emit: metrics
    path "mikado.loci.scores.tsv", emit: scores
    path "*pick.log", emit: log

  script:
    def protein  =   params.protein  ? "-bt ${params.protein}" : ""
    def prot     =   params.protein ? "--xml mikado_prepared.blast.tsv --blast_targets ${params.protein}" : ""
    def junc     =   params.short_reads  ? "--junction ${junction}": ""

    """
    mikado pick --json-conf configuration.yaml --subloci-out mikado.subloci
    """
}
