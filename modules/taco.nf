process taco {

  label "small_job"

  conda "$baseDir/conda-envs/taco-env.yaml"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*taco.gtf'

  input:
    path gtf
    val name

  output:
    path "*taco.gtf", emit: gtf

  script:

    """
    ls -1 *.gtf > gtf_files.txt

    taco_run \
    --num-processes $task.cpus \
    -o taco \
    --filter-min-expr 0 \
    gtf_files.txt

    cp taco/assembly.gtf ./${name}_taco.gtf
    """
}
