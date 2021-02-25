process tama {

  label "long_job"

  conda "$baseDir/conda-envs/tama-env.yaml"
  //container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: '*.gtf'

  input:
    path genome
    path bam

  output:
    path "*.gtf", emit: gtf

  script:

    """
    git clone https://github.com/GenomeRIK/tama

    python ./tama/tama_collapse.py \
    -s $bam -b BAM \
    -f $genome \
    -p ${bam.baseName}_tama

    python ./tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
    ${bam.baseName}_tama.bed \
    ${bam.baseName}_tama.gtf
    """
}
