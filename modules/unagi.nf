process unagi {
  tag "$sample_id"
  label "long_job"

  conda "$baseDir/conda-envs/unagi-env.yaml"
  // container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/transcript_predictions/", mode: 'copy', pattern: 'unagi/*.bed'

  input:
    path genome
    tuple val(sample_id), file(reads)

  output:
    path "unagi/Splicing_Isoforms.bed", emit: gtf

  script:
    def stranded  = params.ont_stranded        ? "--stranded"                   : ""
    def unzip     = reads.Extension == 'gz'    ? "gunzip -f -d -q $reads"       : ""
    def file      = reads.Extension == 'gz'    ? "$reads.baseName"              : "$reads"
    """

    wget https://github.com/iMetOsaka/UNAGI/archive/refs/heads/master.zip
    unzip master.zip
    chmod a+x UNAGI-master/unagi

    $unzip

    cd UNAGI-master

    unagi \
    --input ../$file \
    --genome ../$genome \
    --output ../unagi \
    $stranded

    """
}
