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
    """
    extension=`echo $file | awk -F"." '{ print $NF }'`

    if [ $extension == 'bam' ]
    then
      agat_convert_bed2gff.pl --bed $file -o ${file.getSimpleName}_agat.gff3
    else
      agat_convert_gxf2gxf.pl --gff $file -o ${file.getSimpleName}_agat.gff3
    fi
    """
}
