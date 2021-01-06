process hypo {
  label "assembly"

  conda "$baseDir/conda-envs/hypo-env.yaml"

  publishDir "${params.outdir}/short_polished", mode: 'copy'

  input:
    path genome
    tuple val(sample_id), file(short_reads)
    path long_reads
    val genome_size

  output:
    path "*hypo.fasta", emit: assembly

  script:

    """
    # short read
    minimap2 -t ${task.cpus} -ax sr $genome ${short_reads[1]} ${short_reads[2]} | samtools sort -@ $task.cpus -O BAM- > short.bam

    #long (ONT)
    minimap2 -t $task.cpus -ax map-ont $genome nanopore.fq.gz | samtools sort -@ \$task.cpus -O BAM- > long.bam

    avg_depth = samtools depth short.bam  |  awk '{sum+=\$3} END { print sum/NR}'

    echo -e "${short_reads[1]}\n${short_reads[2]}" > names.txt

    hypo \
    --draft $genome \
    --reads-short names.txt \
    --size-ref $genome_size \
    --coverage-short /$avg_depth \
    --processing-size 96 \
    --bam-sr short.bam \
    --bam-lr long.bam \
    --threads ${task.cpus} \
    --output ${genome.simpleName}_hypo.fasta
    """
}
