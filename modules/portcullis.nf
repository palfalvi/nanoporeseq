process portcullis {

  label "small_job"

  conda "$baseDir/conda-envs/annotation-env.yaml"

  publishDir "${params.outdir}/junction_prediction/", mode: 'copy', pattern: '*_portcullis/junctions.bed'

  input:
    path genome
    tuple file(bam), file(baidx)

  output:
    path "*_portcullis/junctions.bed", emit: junctions

  script:

    def strand = params.orientation ? "--orientation FR" : ""

    """
    portcullis full \
    -t $task.cpus \
    -o ${bam.simpleName}_portcullis \
    $strand \
    $genome \
    $bam
    """
}
