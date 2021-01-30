process hisat2_align {
	tag "$sample_id"
	label 'small_plus'
	publishDir "${params.outdir}/hisat2", mode: 'copy'

  conda "$baseDir/conda-envs/hisat-env.yaml"

  input:
    path genome_idx
    tuple val(sample_id), file(reads)
  output:
    path "*.hisat.bam", emit: bam

  script:
		def strandedness = params.orientation  ? "--fr" : "--rf"
    """
    hisat2 \
		--dta-cufflinks \
		-p $task.cpus \
		-x $genome_idx \
		-1 ${reads[0]} -2 ${reads[1]} |
		samtools sort -@ $task.cpus -O BAM - > ${sample_id}.hisat.bam


		"""
}
