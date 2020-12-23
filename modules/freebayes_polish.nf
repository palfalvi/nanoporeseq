process freebayes {

  label 'assembly'

  conda "$baseDir/conda-envs/freebayes-env.yaml"

  publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    path assembly
    tuple val(sample_id), file(reads)

  output:
    path "*.sample.bam"
    path "*.bam.bai"

  script:
    """
    # bwa index
    bwa index $assembly

    # bwa map
    bwa mem -t ${task.cpus} ${assembly} ${reads[1]} ${reads[2]} > ${assembly}.sam

    samtools sort -@ ${task.cpus} -O bam -o ${assembly}.bam -T ${assembly}.tmp ${assembly}.sam && rm ${assembly}.sam

    samtools index ${assembly}.bam

    # Consensus
    # concat_list.txt:
    #/mnt/gpfsB/scratch/peegee/freebayes/output/tig00000008:1-65284.bcf
    #/mnt/gpfsB/scratch/peegee/freebayes/output/tig00000009:1-422257.bcf
    #/mnt/gpfsB/scratch/peegee/freebayes/output/tig00000011:1-146984.bcf
    #/mnt/gpfsB/scratch/peegee/freebayes/output/tig00000013:1-142623.bcf
    #/mnt/gpfsB/scratch/peegee/freebayes/output/tig00000014:1-193717.bcf
    #/mnt/gpfsB/scratch/peegee/freebayes/output/tig00000015:1-258742.bcf
    #/mnt/gpfsB/scratch/peegee/freebayes/output/tig00000016:1-181964.bcf
    #/mnt/gpfsB/scratch/peegee/freebayes/output/tig00000017:1-137101.bcf
    # ...

    bcftools concat -nf concat_list.txt | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f $fasta -o ${assembly}.bcf
    bcftools index ${assembly}.bcf

    bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f $fasta ${assembly}.bcf > ${assembly}.fasta

    cp ${assembly}.fasta /mnt/gpfsA/home/peegee/cephalotus/assemblies/canu/3_ilmn_polish/${assembly}.fasta
    """
}
