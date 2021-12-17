process agat_extractor {
  tag "$file"
  label "small_job"

  conda "$baseDir/conda-envs/agat-env.yaml"
  //container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/final_transcripts/", mode: 'copy', pattern: '*.fasta'

  input:
    path file
    path genome

  output:
    path "*rna.fasta", emit: mrna
    path "*cds.fasta", emit: cds
    path "*pep.fasta", emit: pep
    path "*promoter2kb.fasta", emit: promoter


  script:
    """
    agat_sp_extract_sequences.pl -gff $file --fasta $genome --mrna -o ${file.baseName}_rna.fasta
    agat_sp_extract_sequences.pl -gff $file --fasta $genome -t cds -o ${file.baseName}_cds.fasta
    agat_sp_extract_sequences.pl -gff $file --fasta $genome --protein --clean_final_stop -o ${file.baseName}_pep.fasta
    agat_sp_extract_sequences.pl -gff $file --fasta $genome -t gene --upstream 2000 -o ${file.baseName}_promoter2kb.fasta
    """
}
