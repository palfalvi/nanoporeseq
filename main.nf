#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
log.info """
=======================================================
       Nanopore Sequence Manipulation pipeline
         https://github.com/palfalvi/rnaseq
=======================================================

Usage:


Mandatory arguments:
      --mode                         Running mode. Supported: basecalling, assembly, annotation, expression.

Basecalling mode:
      --seq_file                     Sequence Files in a csv document. First column is sample name and second column is path to samples (compressed .tar.gz).
      --flowcell
      --kit
      --config_file

Assembly mode:
      --assembler                    One of the following: canu, masurca, flye, miniasm, shasta, haslr
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

// Reading in files
if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Samplesheet file not specified!" }

if ( params.mode == 'basecalling') {

  def sample = file(params.seq_file)

  Channel
    .fromPath( params.seq_file )
    .splitCsv(header:false)
    .map{ row-> tuple( row.sample_id, file(row.read) ) }
    .set { sample_ch }

  log.info 'Found $sample_ch.countLines() samples. '

  if ( config_file ) {
    // Run guppy with config file
    guppy_basecalling(params.guppy, sample_ch)
  } else {
    // Run guppy with kit and flowcell info
  }
}
else if ( params.mode == 'assembly' ) {

}
else if ( params.mode == 'annotation' ) {

} else if ( params.mode == 'expression' ) {

} else if ( !params.mode ) {
  error 'Running mode is not supplied. Please specify with --mode'
} else {
  error 'Invalid running method: $params.mode. Currently supported modes: basecalling, assembly, annotation, expression.'
}
