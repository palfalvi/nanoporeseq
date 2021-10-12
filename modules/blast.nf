process blast {
  label "small_job"

  conda "$baseDir/conda-envs/blast-env.yaml"

  input:
    path reference
    path query
    val blast

  output:
    path '*.blast_sub.tsv', emit: blast

  script:
    """
    $blast -max_target_seqs 5 -num_threads ${task.cpus} -query $query -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop' -db ${reference.baseName[0]} -evalue 0.000001 -out ${reference.baseName}.blast_sub.tsv 2> blast.log
    """
}
