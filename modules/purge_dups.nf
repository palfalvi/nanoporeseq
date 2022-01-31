process purge_dups {

  label 'long_job'

  container "https://depot.galaxyproject.org/singularity/purge_dups:1.2.5--h5bf99c6_1"

  publishDir "${params.outdir}/purge_dups", mode: 'copy'

  input:
  path genome
  path long_reads
  val platform // map-ont (ONT) map-pb (PB CLR) or asm20 (HiFi)

  output:
    path "purged.fa ", emit: purged
    path "hap.fa ", emit: haplo

  script:
    """
    minimap2 -t $task.cpus -x $platform -d ${genome.simpleName}.idx $genome

    minimap2 -x $platform -t $task.cpus ${genome.simpleName}.idx $long_reads | gzip -c - > ${genome.simpleName}_${long_reads.simpleName}.paf.gz

    minimap2 -xasm5 -DP -t $task.cpus $genome $genome | gzip -c - > ${genome.simpleName}_self.paf.gz

    pbcstat ${genome.simpleName}_${long_reads.simpleName}.paf.gz

    calcuts PB.stat > cutoffs 2>calcults.log

    purge_dups -2 -T cutoffs -c PB.base.cov ${genome.simpleName}_self.paf.gz > dups.bed 2> purge_dups.log

    get_seqs dups.bed $genome > purged.fa 2> hap.fa
    """
}
