process trinity_gg {

  label "small_job"

  conda "$baseDir/conda-envs/annotation-env.yaml"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*.gtf'

  input:
    path genome
    tuple file(bam), file(baidx)

  output:
    path "*trinity.fasta", emit: fasta
    path "*.gtf", emit: gtf

  script:

    def strand = params.orientation ? "--SS_lib_type FR" : ""
    // bam2gtf.py is a funannotate module
    """
    Trinity \
    --genome_guided_bam $bam \
    --min_contig_length 100 \
    --genome_guided_max_intron 10000 \
    --max_memory $task.memory \
    --CPU $task.cpus \
    --output trinity_gg \
    --full_cleanup \
    $orientation

    mv ./trinity_gg/Trinity-GG.fasta ./${bam.simpleName}_trinity.fasta

    minimap2 -ax splice:hq --cs -uf $genome ${bam.simpleName}_trinity.fasta > ${bam.simpleName}_trinity.bam

    bam2gtf.py ${bam.simpleName}_trinity.bam ${bam.simpleName}_trinity.gtf
    """
}
