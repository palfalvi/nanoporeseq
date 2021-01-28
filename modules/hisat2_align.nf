process hisat2_align {
	tag "$sample_id"
	label 'small_plus'
	publishDir "${params.out}/hisat2", mode: 'copy'

  conda "$baseDir/conda-envs/hisat-env.yaml"

  input:
    path genome_idx
    tuple val(sample_id), file(reads)
  output:
    path "*.hisat.bam", emit: bam
		path "*.hisat.bam.bai", emit: baidx
  script:
		def strandedness = params.strandedness  ? "--fr" : "--rf"
    """
    hisat2 \
		--dta-cufflinks \
		-p $task.cpus \
		-x $genome_idx \
		-1 ${reads[1]} -2 ${reads[2]} |
		samtools sort -@ $task.cpus -O BAM - > ${sample_id}.hisat.bam

		samtools index ${sample_id}.hisat.bam

		"""
}
