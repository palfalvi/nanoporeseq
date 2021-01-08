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
      --mode                         Running mode. Supported: basecalling, cleanup, assembly, genome_check, annotation, expression.

Universal arguments
      --outdir                       Output directory name. [results]
      -profile                       Sets the running environment. Default is NIBB-BIAS5 PBSPro. 'cde' and 'local' are available to run on NIBB-CDE server or on local machine.
      -bg                            Run the pipeline in the background.
      -resume                        Resume interrupted run and try to catch previously finished processes.

Basecalling mode:
      --seq_file                     Sequence Files in a csv document. First column is 'sample_id' and second column is 'reads' (compressed .tar.gz).
      --flowcell
      --kit
      --config_file                  Guppy config file

Assembly mode:
      --assembler                    One of the following: canu, masurca, flye, miniasm, shasta, haslr, raven.
      --fastq                        Long read fastq file
      --short                        Short read fastq file. Used for some assemblers (e.g masurca, haslr) and for optional polishing.
      --polish                       True/False or software?
      --genome_size                  Expected size of genome.

Assembler specific options?
      --masurca_file                 Setup file for masurca.
      ...

Annotation mode
      --genome                       Reference genome. If not provided, reference-free annotation is attempted.

      --lorean_long                  Long read fastq file.
      --lorean_short                 Short read fastq file. If pair end, write them seprataed by a come, e.g. "read1.fastq,read2.fastq".
      --lorean_proteins              Protein homologs for the gene prediction.
      --lorean_iproscan              Boolean. If provided, LoReAn runs with InterProScan.
      --lorean_adapters              Adapter fasta file.
      --lorean_species               Species name.
      --lorean_prefix                Prefix for the gene names.

""".stripIndent()
}

params.help = false
if (params.help){
helpMessage()
exit 0
}

// Include basecallers
include { guppy_basecalling } from './modules/guppy_basecalling.nf'

//Include cleanup tools
include { nanolyse } from './modules/nanolyse.nf'
include { pychopper } from './modules/pychopper.nf'

// Include assemblers
include { canu } from './modules/canu.nf'
include { masurca } from './modules/masurca.nf'
include { raven } from './modules/raven.nf'
include { flye } from './modules/flye.nf'
include { miniasm } from './modules/miniasm.nf'
include { nextdenovo } from './modules/nextdenovo.nf'
include { wtdbg } from './modules/wtdbg.nf'


// Include polishing tools
include { minimap2 as minimap2_1; minimap2 as minimap2_2; minimap2 as minimap2_3 } from './modules/minimap.nf'
include { bwa_index } from './modules/bwa_index.nf'
include { bwa_mem } from './modules/bwa_mem.nf'
include { bam_coverage } from './modules/bam_coverage.nf'
include { racon as racon1; racon as racon2; racon as racon3 } from './modules/racon.nf'
include { medaka } from './modules/medaka.nf'
include { hypo } from './modules/hypo.nf'
include { freebayes_bwa } from './modules/freebayes-bwa_map.nf'
include { freebayes_call } from './modules/freebayes-call.nf'
include { freebayes_consensus } from './modules/freebayes-consensus.nf'

// Include annotation tools
include { lorean } from './modules/lorean.nf'

// Include QC tools
include { busco as busco_vir; busco as busco_emb; busco as busco_eud } from './modules/busco.nf'
include { quast } from './modules/quast.nf'
include { multiqc } from './modules/multiqc.nf'
include { kat } from './modules/kat.nf'



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
else if ( params.mode == 'cleanup' ) {
  log.info 'Starting read clean-up'


    reads = Channel.fromPath( params.fastq )
    reads.subscribe {  println "Reads provided: $it"  }

    if ( !params.rna ) {
      nanolyse( reads )
    }
    else {
      pychopper( reads )
    }

}
else if ( params.mode == 'assembly' ) {


  if ( params.assembly && !params.assembler ) {
    // assembly file is provided, but no assembler software.
    // Useful for only polishing or QC.
    log.info "Skipping initial assembly as it is provided: $params.assembly "
    assembly = params.assembly
  }
  else {
    log.info "Starting assembly protocol with $params.assembler ... "
  }


  // Canu, out.assembly, out.gfa
  if ( params.assembler == 'canu' ) {

  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."
  //  params.genome_size ? log.info "Estimated genome size: $it" : error "Estimated genome size is missing but needed for canu. Please provide with --genome_size."

    canu(params.fastq, params.genome_size)

    assembly = canu.out.assembly
    gfa = canu.out.gfa
  }

  // MaSuRCa, out.assembly
  if ( params.assembler == 'masurca' ) {

  //  params.masurca_file ? log.info "MaSuRCa config file: $it" : error "MaSuRCa config file is not provided. Please specify with --masurca_file parameter."

    masurca( params.masurca_file )

    assembly = masurca.out.assembly
  }

  // Flye: out.assembly, out.gfa
  if ( params.assembler == 'flye' ) {

  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."
  //  params.genome_size ? log.info "Estimated genome size: $it" : error "Estimated genome size is missing but needed for flye. Please provide with --genome_size."

    flye(params.fastq, params.genome_size)

    assembly = flye.out.assembly
    gfa = flye.out.gfa
  }

  // Miniasm out.assembly, out.gfa
  if ( params.assembler == 'miniasm' ) {

  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."

    miniasm(params.fastq)

    assembly = miniasm.out.assembly
    gfa = miniasm.out.gfa
  }

  // wtdbg2: out.assembly
  if ( params.assembler == 'wtdbg' || params.assembler == 'wtdbg2' || params.assembler == "redbean" ) {

  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."
  //  params.genome_size ? log.info "Estimated genome size: $it" : error "Estimated genome size is missing but needed for wtdbg2. Please provide with --genome_size."

    wtdbg2(params.fastq, params.genome_size)

    assembly = wtdbg2.out.assembly
  }

  // NextDenovo: out.assembly
  if ( params.assembler == 'nextdenovo') {

  //  params.nextdenovo_bin ? log.info "Using NextDenovo: $it" : error "NextDenovo executable file is not provided. Please note that NextDenovo cannot be installed in a conda environment, thus you need to provide the local executable location."
  //  params.nextdenovo_cfg ? log.info "Using NextDenovo config file: $it" : error "NextDenovo config file is not provided. Please provide with --nextdenovo_cfg"

    nextdenovo(params.nextdenovo_bin, params.nextdenovo_cfg)

    assembly = nextdenovo.out.assembly
    gfa = nextdenovo.out.gfa
  }

  // RAVEN assembler: out.assembly out.gfa
  if ( params.assembler == 'raven') {

    // params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."

    raven(params.fastq)

    assembly = raven.out.assembly
    gfa = raven.out.gfa

  }


  if ( params.polish ) {

    minimap2_1(params.fastq, assembly)
    racon1(params.fastq, minimap2_1.out.map, assembly)

    medaka(params.fastq, racon1.out.assembly)

    //medaka.out.assembly

    assembly = medaka.out.assembly

  }

  if ( params.short_polish ) {

    if ( !params.short_reads ) {
      error 'Short reads are not provided. Please provide short reads as --short_reads /path/to/short.fastq'
    } else {
      // Short read mapping

      short_r = Channel.fromFilePairs( params.short_reads )

      short_r.subscribe {  println "Short reads provided: $it"  }

      if ( params.short_polish_map == "bwa" ) {
        // Mapping with bwa-mem

        bwa_index( assembly )

        bwa_mem( short_r, assembly, bwa_index.out.index )

        short_bam = bwa_mem.out.bam
        short_baidx = bwa_mem.out.baidx

        bam_coverage( short_bam )

        bam_coverage.out.coverage.subscribe( println "Short read coverage is $it x." )

      }
      else if ( params.short_polish_map === "minimap2" ) {
        // Mapping with minimap2 -ax sr
      }
      else {
        error 'Unknown mapping method: ${params.short_polish_map}. Please choose from bwa, minimap2 or contact developers.'
      }

    }

    if ( params.short_polish == 'freebayes' | params.short_polish == 'vgp' | params.short_polish == true ) {
      // Verterae Genome Project polishing pipeline with freebayes and bcftools

      // Split genome by contigs
      //Channel.fromPath( assembly )
      //  .splitFasta( by: 100, file: true )
      //  .set { contigs_ch }
      // Run freebayes on each contigs

      contig_index = Channel.from(1..100)

      freebayes_call( assembly, bam_coverage.out.coverage, short_bam, short_baidx, contig_index )

      freebayes_call.out.collectFile(name: 'concat_list.txt', newLine: true, sort: true).set { bcf_list }

      freebayes_consensus( assembly, bcf_list )

      polished_assembly = freebayes_consensus.out.assembly

    }
    else if ( params.short_polish == 'nextpolish' ) {
      // Nextpolish polishing
      if ( !params.nextpolish ) {
        error 'NextPolish source is not provided. Please install NextPolish locally and provide as --nextpolish /PATH/TO/nextpolish or consider using another polishing method (e.g. vgp).'
      }

    }
    else if ( params.short_polish == 'pilon' ) {
      // pilon polishing


      pilon( assembly, short_bam, short_baidx)

      polished_assembly = pilon.out.assembly
    }
    else if ( params.short_polish == 'polca' ) {
      // POLCA polishing
    }
    else if ( params.short_polish == 'hypo' ) {
      // HyPo polishing

      hypo( assembly, short_r, params.fastq, params.genome_size )
      polished_assembly = hypo.out.assembly
    }

  }
  else if ( params.polish ) {
    polished_assembly = medaka.out.assembly
  } else {
    polished_assembly = assembly
  }



  if ( !params.skip_qc) {
    // QC
    quast(polished_assembly)
    busco_eud(polished_assembly, "eudicots_odb10", "genome")
    busco_emb(polished_assembly, "embryophyta_odb10", "genome")
    busco_vir(polished_assembly, "viridiplantae_odb10", "genome")

    multiqc(quast.out.summary.mix(busco_eud.out, busco_emb.out, busco_vir.out).collect(), "$baseDir/${params.outdir}")

    if ( params.short_reads ) {

      short_reads = Channel.fromFilePairs( params.short_reads )

      kat(polished_assembly, short_reads)
    }
  }



// deduplication?




}
else if ( params.mode == 'genome_check' ) {

  // Run quast and busco on an assembled genome
  quast(params.genome)
  busco_eud(params.genome, "eudicots_odb10", "genome")
  busco_emb(params.genome, "embryophyta_odb10", "genome")
  busco_vir(params.genome, "viridiplantae_odb10", "genome")

  multiqc(quast.out.summary.mix(busco_eud.out, busco_emb.out, busco_vir.out).collect(), "$baseDir/${params.outdir}")

  if ( params.short_reads ) {

    short_reads = Channel.fromFilePairs( params.short_reads )

    kat(polished_assembly, short_reads)
  }

}

else if ( params.mode == 'annotation' ) {
  log.info "Starting annotation protocol ... "

  if ( params.genome != false ) {
    // Genome file is provided, run LoReAn
    log.info "Genome file provided: ${params.genome}"
    log.info "LoReAn annotation is starting ..."

    lorean(params.genome, params.lorean_proteins)

  } else {
    log.info 'No reference genome is provided for transcript annotation.'
    log.info 'Attempting de novo transcript assembly...'
    // No genome file is provided, do de novo transcript assembly
    log.info 'Sorry, these functions are not yet implemented .... '
    if ( !params.short_rna && !params.long_rna ) {
      // Both short and long reads are provided.
      //rnaSPAdes?
    }
    else if( !params.long_rna ) {
      // Only long RNA-seq data provided
      // IsOnClust2?
    }
    else if ( !params.short_rna ) {
      //Only short RNA is provided
      //Trinity?
      //rnaSPAdes?
    }
    else { error 'No read file is provided for transcriptome assembly. Please provide fastq files with --short_rna and/or --long_rna options.'}
  }



}
else if ( params.mode == 'expression' ) {
  log.info "Starting gene expression protocol ... "
  log.info 'Except... I cannot do that yet. Sorry.'
}
else if ( !params.mode ) {
  error 'Running mode is not supplied. Please specify with --mode'
} else {
  error 'Invalid running method: $params.mode. Currently supported modes: basecalling, assembly, annotation, expression.'
}


}
