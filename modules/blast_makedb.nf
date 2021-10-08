process blast_makedb {
  label "small_job"

  conda "$baseDir/conda-envs/blast-env.yaml"

  input:
    path reference
    value dbtype

  output:
    path('*')

  script:
    def unzip      = reference.Extension == 'gz'    ? "gunzip -f -d -q $reference"       : ""
    def file       = reference.Extension == 'gz'    ? "$reference.baseName"              : "$reference"
    """
    $unzip

    makeblastdb -in $file -dbtype $dbtype -parse_seqids > blast_prepare.log
    """
}
