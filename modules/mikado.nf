process mikado {
  label "long_job"

  conda "$baseDir/conda-envs/mikado-env.yaml"
  // container "peegee/nanoporeseq:latest"

  publishDir "${params.outdir}/mikado/", mode: 'copy', pattern: 'mikado*'

  input:
    path genome
    path('*')
    // path gtf_file
    path scoring
    path junction

  output:
    path "mikado_prepared.gtf", emit: gtf
    path "mikado_prepared.fasta", emit: fasta
    path "prepare.log", emit: log

  script:
    def sh  = params.short_reads    ? "${projectDir}/scripts/short_gtf.txt"   : ""
    def ont = params.ont_reads      ? "${projectDir}/scripts/ont_gtf.txt"     : ""
    def pb  = params.pb_reads       ? "${projectDir}/scripts/pb_gtf.txt"      : ""
    def pr  = !params.skip_abinitio ? "${projectDir}/scripts/prot_gtf.txt"    : ""

    def protein  =   params.protein  ? "-bt ${params.protein}" : ""
    def blastdb  =   params.protein ? "makeblastdb -in ${params.protein} -dbtype prot -parse_seqids > blast_prepare.log" : ""
    def blastjob =   params.protein ? "blastx -max_target_seqs 5 -num_threads 10 -query mikado_prepared.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop' -db ${params.protein} -evalue 0.000001 -out mikado_prepared.blast.tsv 2> blast.log " : ""
    def prot     =   params.protein ? "--xml mikado_prepared.blast.tsv" : ""
    def junc     =   junction != []  ? "--junction ${junction}": ""

    """
    cat $sh $ont $pb $pr > file.txt
    sed -e 's/ /\t/g' file.txt | sed 's/.*/${params.outdir}&/' > gtf_list.txt

    mikado configure \
    --list gtf_list.txt \
    --reference $genome \
    --mode permissive \
    --scoring $scoring  \
    $junc \
    $protein \
    --threads $task.cpus \
    configuration.yaml

    mikado prepare \
    --json-conf configuration.yaml

    $blastdb
    $blastjob

    #TransDecoder

    mikado serialise --json-conf configuration.yaml $prot --orfs mikado.bed --blast_targets
    mikado pick --json-conf configuration.yaml --subloci_out mikado.subloci.gff3
    """
}
