process agat_converter {

  label "small_job"

  conda "$baseDir/conda-envs/agat-env.yaml"
  //container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/agat_transcripts/", mode: 'copy', pattern: '*_agat.gff3'

  input:
    path file

  output:
    path "*_agat.gff3", emit: gff

  script:
    def agat = file.Extension == 'bed' ? "agat_convert_bed2gff.pl --bed $file" : "agat_convert_sp_gxf2gxf.pl --gff $file"

    """
    $agat -o ${file.SimpleName}_agat.gff3
    """
}
