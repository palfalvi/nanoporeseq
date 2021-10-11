process mikado_prepare {
  label "long_job"

  conda "$baseDir/conda-envs/mikado-env.yaml"
  // container "peegee/nanoporeseq:latest"

  // publishDir "${params.outdir}/mikado/", mode: 'copy', pattern: 'mikado*'

  input:
    path genome
    path('*')
    path scoring
    //path junction

  output:
    path "configuration.yaml", emit: config
    path "mikado_prepared.fasta", emit: fasta
    path "mikado_prepared.gtf", emit: gtf

  script:
    def sh  = params.short_reads    ? "${projectDir}/scripts/short_gtf.txt"   : ""
    def ont = params.ont_reads      ? "${projectDir}/scripts/ont_gtf.txt"     : ""
    def pb  = params.pb_reads       ? "${projectDir}/scripts/pb_gtf.txt"      : ""
    def pr  = !params.skip_abinitio ? "${projectDir}/scripts/prot_gtf.txt"    : ""

    def protein  =   params.protein  ? "-bt ${params.protein}" : ""
    def junc     =   params.short_reads  ? "--junction ${junction}": ""
    """
    cat $sh $ont $pb $pr > file.txt
    sed -e 's/ /\t/g' file.txt > gtf_list.txt

    mikado configure \
    --list gtf_list.txt \
    --reference $genome \
    --mode permissive \
    --scoring $scoring  \
    $junc \
    $protein \
    --threads $task.cpus \
    configuration.yaml

    mikado prepare \
    --json-conf configuration.yaml
    """
}
