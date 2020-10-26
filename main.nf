#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
=======================================================
       Nanopore Sequence Manipulation pipeline
         https://github.com/palfalvi/rnaseq
=======================================================

Usage:


Mandatory arguments:
      --mode                         Running mode. Supported: basecalling, assembly, annotation, expression.

Basecalling mode:
      --seq_fofn                     Sequence File of File Names. One compressed file per line.
      --flowcell
      --kit
      --barcode_kit

Assembly mode:
      --assembler                    One of the following: canu, masurca, flye, miniasm, haslr
      --fastq                        Long read fastq file
      --short                        Short read fastq file. Used for some assemblers (e.g masurca, haslr) and for optional polishing.
      --polish                       True/False or software
      --genome_size                  Expected size of genome.
      --busco                        BUSCO species (list?)

Assembler specific options?
      --masurca_file                 Setup file for masurca.
      ...

Annotation mode
      --fastq                        Long read fastq file.
      --genome                       Reference genome. If not provided, reference-free annotation is attempted.
      --busco                        BUSCO species (lsit?)
      --

Computer allocation settings
      -profile                       Sets the running environment. Default is NIBB-BIAS5 PBSPro. 'cde' and 'local' are available to run on NIBB-CDE server or on local machine.
""".stripIndent()
}

params.help = false
if (params.help){
helpMessage()
exit 0
}


if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Samplesheet file not specified!" }
