process edta_softmask {

  label "long_job"

  conda "$baseDir/conda-envs/edta-env.yaml"

  publishDir "${params.outdir}/soft_masking", mode: 'copy'

  input:
    path genome

  output:
    path "*.mod.EDTA.TElib.fa", emit: novel_tes
    path "*.mod.EDTA.TEanno.gff", emit: te_anno
    path "*.mod.EDTA.TEanno.sum", emit: summary
    path "*.mod.MAKER.masked", emit: masked // This can be used for masking

  script:
    def cds       = params.cds  ? "--cds ${params.cds}"        : ""

    """
    EDTA.pl \
    --genome $genome \
    $cds \
    --sensitive 1 \
    --anno 1 \
    --threads $task.cpus
    """
}
