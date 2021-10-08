process mikado_serialise {
  label "long_job"

  conda "$baseDir/conda-envs/mikado-env.yaml"
  // container "peegee/nanoporeseq:latest"

  // publishDir "${params.outdir}/mikado/", mode: 'copy', pattern: 'mikado*'

  input:
    path genome
    path config
    path blastp
    path transdecoder
    path scoring
    path junction

  output:
    path "mikado.db", emit: db

  script:
    def prot     =   params.protein ? "--xml mikado_prepared.blast.tsv --blast_targets ${params.protein}" : ""
    def junc     =   params.short_reads  ? "--junction ${junction}": ""

    """
    mikado serialise --json-conf $config $prot --orfs $transdecoder $junc
    """
}
