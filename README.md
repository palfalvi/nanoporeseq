# nanoporeseq

  THIS PROJECT IS EXPERIMENTAL!

nextflow pull palfalvi/nanoporeseq


nextflow run palfalvi/nanoporeseq -r main --fastq '/mnt/gpfsA/home/peegee/cephalotus/ont_genome/191209_cep_ont_reads.fastq' --mode assembly --assembler raven --polish --outdir raven_assembly

nextflow run palfalvi/nanoporeseq -r main --fastq '/mnt/gpfsA/home/peegee/cephalotus/ont_genome/191209_cep_ont_reads.fastq' --mode assembly --assembler miniasm --polish --outdir miniasm_assembly

nextflow run palfalvi/nanoporeseq -r main --fastq '/mnt/gpfsA/home/peegee/cephalotus/ont_genome/191209_cep_ont_reads.fastq' --mode assembly --assembler canu --polish --genome_size 2g --outdir canu_assembly

nextflow run palfalvi/nanoporeseq -r main --fastq '/mnt/gpfsA/home/peegee/cephalotus/ont_genome/191209_cep_ont_reads.fastq' --mode assembly --assembler flye --polish --genome_size 2g --outdir flye_assembly

nextflow run palfalvi/nanoporeseq -r main --fastq '/mnt/gpfsA/home/peegee/cephalotus/ont_genome/191209_cep_ont_reads.fastq' --mode assembly --assembler wtdbg2 --polish --genome_size 2g --outdir wtdbg2_assembly


// genome check
nextflow run palfalvi/nanoporeseq -r main --mode genome_check --genome /mnt/gpfsB/scratch/peegee/nextflow/nanopore/assembly/raven_assembly/raven_assembly.fasta --outdir qc

GeneMarkES needs to be installed for annotation !!!
need to do:
cd ~/bin
change_path_in_perl_scripts.pl "/usr/bin/env perl"
Otherwise braker will not look into conda perl modules
