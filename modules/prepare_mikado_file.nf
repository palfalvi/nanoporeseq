process prepare_mikado_file {

  label "small_job"

  publishDir "${params.outdir}/mikado/", mode: 'copy', pattern: 'gtf_list.txt'

  input:
    path gtf

  output:
    path "gtf_list.txt", emit: gtf_list

  script:
    def sh  = params.short_reads    ? "${projectDir}/scripts/short_gtf.txt"   : ""
    def ont = params.ont_reads      ? "${projectDir}/scripts/ont_gtf.txt"     : ""
    def pb  = params.pb_reads       ? "${projectDir}/scripts/pb_gtf.txt"      : ""
    def pr  = !params.skip_abinitio ? "${projectDir}/scripts/prot_gtf.txt"    : ""
    """
    cat $sh $ont $pb $pr > file.txt
    sed -e 's/ /\t/g' file.txt > gtf_list.txt 
    """
}
