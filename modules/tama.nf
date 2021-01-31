process tama {

  label "small_plus"

  conda "$baseDir/conda-envs/tama-env.yaml"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*.gtf'

  input:
    path genome
    tuple file(bam), file(baidx)

  output:
    path "*.gtf", emit: gtf

  script:

    """
    git clone https://github.com/GenomeRIK/tama

    python ./tama/tama_collapse.py \
    -s $bam -b BAM \
    -f $genome \
    -p ${bam.baseName}_tama \
    -sj sj_priority \
    -lde 5 -sjt 20 -a 100 -z 100

    python ./tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
    ${bam.baseName}_tama.bed \
    ${bam.baseName}_tama.gtf
    """
}
