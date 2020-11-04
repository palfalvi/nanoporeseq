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
      --seq_file                     Sequence Files in a csv document. First column is 'sample_id' and second column is 'reads' (compressed .tar.gz).
      --flowcell
      --kit
      --config_file                  Guppy config file

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
      --busco                        BUSCO species (list?)
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

// Include basecallers
include { guppy_basecalling } from './modules/guppy_basecalling.nf'

// Include assemblers
include { canu } from './modules/canu.nf'
include { masurca } from './modules/masurca.nf'
include { raven } from './modules/raven.nf'
include { flye } from './modules/flye.nf'
include { miniasm } from './modules/miniasm.nf'
include { nextdenovo } from './modules/nextdenovo.nf'
include { wtdbg } from './modules/wtdbg.nf'


// Include polishing tools
include { racon as racon1; racon as racon2; racon as racon3 } from './modules/racon.nf'

// Include QC tools
include { busco as busco_vir; busco as busco_emb; busco as busco_eud } from './modules/busco.nf'
include { quast } from './modules/quast.nf'
include { multiqc } from './modules/multiqc.nf'

workflow {

if ( params.mode == 'basecalling') {

  log.info "Starting basecalling protocol ... "

  def sample = file(params.seq_file)

  Channel
    .fromPath( params.seq_file )
    .splitCsv(header:true)
    .map{ row-> tuple( row.sample_id, file(row.read) ) }
    .set { sample_ch }

  sample_ch.subscribe {  println "ONT library input provided: $it"  }

  //log.info "Found $sample.countLines() samples."

  guppy_basecalling(sample_ch, params.guppy)

}
else if ( params.mode == 'assembly' ) {
  log.info "Starting assembly protocol with $params.assembler ... "

  // Canu, out.assembly, out.gfa
  if ( params.assembler == 'canu' ) {

  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."
  //  params.genome_size ? log.info "Estimated genome size: $it" : error "Estimated genome size is missing but needed for canu. Please provide with --genome_size."

    canu(params.fastq, params.genome_size)

  }

  // MaSuRCa, out.assembly
  if ( params.assembler == 'masurca' ) {

  //  params.masurca_file ? log.info "MaSuRCa config file: $it" : error "MaSuRCa config file is not provided. Please specify with --masurca_file parameter."

    masurca( params.masurca_file )

  }

  // Flye: out.assembly, out.gfa
  if ( params.assembler == 'flye' ) {

  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."
  //  params.genome_size ? log.info "Estimated genome size: $it" : error "Estimated genome size is missing but needed for flye. Please provide with --genome_size."

    flye(params.fastq, params.genome_size)

  }

  // Miniasm out.assembly, out.gfa
  if ( params.assembler == 'miniasm' ) {

  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."

    miniasm(params.fastq)

  }

  // wtdbg2: out.assembly
  if ( params.assembler == 'wtdbg' || params.assembler == 'wtdbg2' || params.assembler == "redbean" ) {

  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."
  //  params.genome_size ? log.info "Estimated genome size: $it" : error "Estimated genome size is missing but needed for wtdbg2. Please provide with --genome_size."

    flye(params.fastq, params.genome_size)

  }

  // NextDenovo: out.assembly
  if ( params.assembler == 'nextdenovo') {

  //  params.nextdenovo_bin ? log.info "Using NextDenovo: $it" : error "NextDenovo executable file is not provided. Please note that NextDenovo cannot be installed in a conda environment, thus you need to provide the local executable location."
  //  params.nextdenovo_cfg ? log.info "Using NextDenovo config file: $it" : error "NextDenovo config file is not provided. Please provide with --nextdenovo_cfg"

    nextdenovo(params.nextdenovo_bin, params.nextdenovo_cfg)

  }

  // RAVEN assembler: out.assemblym out.gfa
  if ( params.assembler == 'raven') {

    // params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."

    raven(params.fastq)

    quast(raven.out.assembly)
    busco_eud(raven.out.assembly, "eudicots_odb10", "genome")
    busco_emb(raven.out.assembly, "embryophyta_odb10", "genome")
    busco_vir(raven.out.assembly, "viridiplantae_odb10", "genome")

    multiqc(quast.out.concat(busco_eud.out, busco_emb.out, busco_vir.out), params.outdir)

  }

 // Polishing?

// deduplication?

// k-mer analysis (KAT)

 // quast

 // busco



}
else if ( params.mode == 'genome_check' ) {

  // Run quast and busco on an assembled genome
  quast(params.genome)
  busco_eud(params.genome, "eudicots_odb10", "genome")
  busco_emb(params.genome, "embryophyta_odb10", "genome")
  busco_vir(params.genome, "viridiplantae_odb10", "genome")

  multiqc(quast.out.concat(busco_eud.out, busco_emb.out, busco_vir.out), params.outdir)
}
else if ( params.mode == 'annotation' ) {
  log.info "Starting annotation protocol ... "
}
else if ( params.mode == 'expression' ) {
  log.info "Starting gene expression protocol ... "
}
else if ( !params.mode ) {
  error 'Running mode is not supplied. Please specify with --mode'
} else {
  error 'Invalid running method: $params.mode. Currently supported modes: basecalling, assembly, annotation, expression.'
}


}
