process minimap2 {

  cpus 32
  memory '1.5T'
  queue 'MPP,smpl'

  conda "$baseDir/conda-envs/minimap-env.yaml"

  // publishDir "${params.outdir}/bwa", mode: 'copy'

  input:
    path fastq
    path assembly

  output:
    path "*.paf", emit: map

  script:
    """
    minimap2 -x map-ont -t ${task.cpus} $assembly $fastq > ${assembly.simpleName}-${fastq.simpleName}.paf
    """
}
