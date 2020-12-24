process hypo {
  label "assembly"

  conda "$baseDir/conda-envs/kat-env.yaml"

  publishDir "${params.outdir}/short_polished", mode: 'copy'

  input:
    path genome
    tuple val(sample_id), file(reads)
    path long_reads
    value genome_size
    value coverage // ? can we calculate somehow?

  output:
    path "*hypo.fasta", emit: assembly

  script:
    def reference  = params.cu_reference        ? "--reference ${params.cu_reference}"  : "--reference $projectDir/conf/ont_control_dna.fasta"
    def qual       = params.cu_qual             ? "--quality ${params.cu_qual}"         : ""
    def length     = params.cu_length           ? "--length ${params.cu_length}"        : ""
    def headcrop   = params.cu_headcrop         ? "--headcrop ${params.cu_headcrop}"    : ""
    def decompress = params.gz                  ? "gunzip -c $reads"                    : "cat $reads"

    """
    # short read
    minimap2 -t ${task.cpus} -ax sr $genome ${short_reads[1]} ${short_reads[2]} | samtools sort -@ $task.cpus -O BAM- > short.bam

    #long (ONT)
    minimap2 -t $task.cpus -ax map-ont $genome nanopore.fq.gz | samtools sort -@ \$task.cpus -O BAM- > long.bam

    avg_depth = samtools depth short.bam  |  awk '{sum+=\$3} END { print sum/NR}' # Calculate average coverage of reads

    echo -e "${short_reads[1]}\n${short_reads[2]}" > names.txt

    hypo \
    --draft $genome \
    --reads-short @names.txt \
    --size-ref $genome_size \
    --coverage-short $avg_depth \
    --processing-size 96 \
    --bam-sr short.bam \
    --bam-lr long.bam \
    --threads ${task.cpus} \
    --output ${genome.simpleName}_hypo.fasta
    """
}
