process agat_longest {
  tag "$file"
  label "small_job"

  conda "$baseDir/conda-envs/agat-env.yaml"
  //container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/final_transcripts/", mode: 'copy', pattern: '*_longestIsoforms.gff3'

  input:
    path file

  output:
    path "_longestIsoforms.gff3", emit: gff

  script:
    """
    agat_sp_keep_longest_isoform.pl -gff $file -o ${file.SimpleName}_longestIsoforms.gff3
    """
}
