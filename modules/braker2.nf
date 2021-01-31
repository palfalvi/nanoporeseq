process braker2 {

  label "long_job"

  conda "$baseDir/conda-envs/braker2-env.yaml"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy'

  input:
    path genome
    path bam

  output:
    path "*braker.gtf", emit: gtf

  script:
    def protein   = params.protein                          ? "--prot_seq=${params.proteins}" : ""
    def mapping   = bam!=null                               ? "--bam $bam"                    : ""
    def tmark     = mapping != ""                           ? "t"                            : ""
    def pmark     = params.protein                          ? "p"                             : ""
    def tpsmark   = (tmark+pmark).length()>0                 ? tmark+pmark                     : "s"
    def species   = params.species                          ?: "sp1"

    """
    braker.pl \
    --genome=$genome \
    $mapping \
    $protein \
    --prg=gth \
    --e${tpsmark}mode \
    --gth2traingenes \
    $species \
    --softmasking \
    --cores $task.cpus
    """
}
