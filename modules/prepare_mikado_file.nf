process prepare_mikado_file {

  label "small_job"

  conda "$baseDir/conda-envs/annotation-env.yaml"

  publishDir "${params.outdir}/mikado/", mode: 'copy', pattern: 'mikado*'

  input:
    path gtf

  output:
    path "*.txt", emit: gtf_list

  script:
    def sh  = params.short_reads    ? "${projectDir}/scripts/short_gtf.txt"   : ""
    def ont = params.ont_reads      ? "${projectDir}/scripts/ont_gtf.txt"     : ""
    def pb  = params.pb_reads       ? "${projectDir}/scripts/pb_gtf.txt"      : ""
    def pr  = !params.skip_abinitio ? "${projectDir}/scripts/prot_gtf.txt"    : ""
    """
    cat $sh $ont $pb $pr > gtf_list.txt
    """
}
