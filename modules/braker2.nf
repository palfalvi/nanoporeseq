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
    if (mapping != "" && params.protein) {
      def mark = "--etpmode"
    } else if (params.protein) {
      def mark = "--epmode"
    } else if (mapping != "") {
      def mark = ""
    } else {
      def mark = "--esmode"
    }
    def species   = params.species                          ?: "sp1"

    """
    braker.pl \
    --genome=$genome \
    $mapping \
    $protein \
    --prg=gth \
    $mark \
    --gth2traingenes \
    $species \
    --softmasking \
    --cores $task.cpus
    """
}
