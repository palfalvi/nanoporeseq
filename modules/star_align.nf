process star_align {
	tag "$sample_id"
	label 'small_plus'
	publishDir "${params.out}/star", mode: 'copy'

  conda "$baseDir/conda-envs/star-env.yaml"

  input:
    path genome_idx
    tuple val(sample_id), file(reads)
  output:
    path "*.bam", emit: bam
		
  script:
    """
    mkdir $sample_id
    STAR \
		 --twopassMode Basic \
		--runThreadN $task.cpus \
		--genomeDir $genome_idx \
		--readFilesIn $reads \
		--readFilesCommand zcat \
		--outFileNamePrefix ${sample_id}.star \
		--outSAMattributes XS \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMstrandField intronMotif
		"""
}
