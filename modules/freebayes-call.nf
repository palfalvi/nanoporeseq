process freebayes_call {
  tag "$contig"
  label 'large_mem'

  conda "$baseDir/conda-envs/freebayes-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    path contig
    val avg_depth
    path bam
    path baidx

  output:
    path "*.bcf"

  script:
    """
    samtools faidx ${contig}

    contig=`sed -n 1p ${contig}.fai | awk '{print \$1}'`
    contig_no_pipe=`echo \$contig | sed 's/|/_/g'`
    end=`sed -n 1p ${contig}.fai | awk '{print \$2}'`

    freebayes --bam $bam --region=\$contig:1-\$end --skip-coverage \$((${avg_depth}*12)) -f ${contig} | bcftools view --no-version -Ou > \${contig_no_pipe}.bcf
    """
}
