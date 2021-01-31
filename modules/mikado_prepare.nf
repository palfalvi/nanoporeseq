process prepare_mikado_file {

  label "small_job"

  conda "$baseDir/conda-envs/annotation-env.yaml"

  publishDir "${params.outdir}/mikado/", mode: 'copy', pattern: 'mikado*'

  input:
    path genome
    path gtf
    path gtf_file
    path scoring
    path junction

  output:
    path "mikado_prepared.gtf", emit: gtf
    path "mikado_prepared.fasta", emit: fasta
    path "prepare.log", emit: log
    path "configuration.yaml", emit: config

  script:
    def protein  =   params.protein  ? "-bt $params.protein" : ""
    def blastdb  =   params.protein ? "makeblastdb -in $params.protein -dbtype prot -parse_seqids > blast_prepare.log" : ""
    def blastjob =   params.protein ? "blastx -max_target_seqs 5 -num_threads 10 -query mikado_prepared.fasta -outfmt 5 -db $params.protein -evalue 0.000001 2> blast.log | sed '/^$/d' | gzip -c - > mikado.blast.xml.gz"
    def prot     =   params.protein ? "--xml mikado.blast.xml.gz" : ""
    def junc     =   junction != ''  ? "--junction $junction": ""

    """
    mikado configure \
    --list $gtf_file \
    --reference $genome \
    --mode permissive \
    --scoring $scoring  \
    $junc \
    $protein \
    --threads = $task.cpus \
    configuration.yaml

    mikado prepare \
    --json-conf configuration.yaml

    $blastdb
    $blastjob

    mikado serialise --json-conf configuration.yaml $prot --orfs mikado.bed --blast_targets
    mikado pick --json-conf configuration.yaml --subloci_out mikado.subloci.gff3
    """
}
