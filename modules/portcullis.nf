process portcullis {

  label "small_job"

  conda "$baseDir/conda-envs/portcullis-env.yaml"

  publishDir "${params.outdir}/junction_prediction/", mode: 'copy', pattern: '*_portcullis/junctions.bed'

  input:
    path genome
    path bam

  output:
    path "*_portcullis/3-filt/portcullis_filtered.pass.junctions.bed", emit: junctions

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
