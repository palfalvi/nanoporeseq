process braker2 {

  label "small_job"

  conda "$baseDir/conda-envs/braker2-env.yaml"

  publishDir "${params.outdir}/braker2/", mode: 'copy'

  input:
    path genome
    path bam

  output:
    path "braker.gtf", emit: gtf

  script:
    def protein   = params.protein                          ? "--prot_seq=${params.proteins}" : ""
    def mapping   = bam!=null                               ? "--bam $bam"                    : ""
    //def shorts    = params.short_reads.or(params.ont_reads) ? "--bam ${bam}"                  : ""
    def tmark     = shorts != ""                            ? "t"                             : ""
    def pmark     = params.protein                          ? "p"                             : ""
    def tpsmark   = (tmark+park).length()>0                 ? tmark+pmark                     : "s"

    """
    braker.pl \
    --genome=$genome \
    $shorts \
    $protein \
    --prg=gth \
    --e${tpsmark}mode \
    --gth2traingenes \
    --species = $params.species \
    --softmasking \
    --cores $task.cpus
    """
}
