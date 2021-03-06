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
      --mode                         Running mode. Supported: basecalling, cleanup, assembly, genome_qc, annotation, expression.

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
include { pilon } from './modules/pilon.nf'

// Include annotation tools
include { edta } from './modules/edta.nf'
include { star_idx } from './modules/star_index.nf'
include { star_align } from './modules/star_align.nf'
include { hisat2_idx } from './modules/hisat2_index.nf'
include { hisat2_align } from './modules/hisat2_align.nf'
include { minimap_rna } from './modules/minimap_rna.nf'
include { bam_merge as merge_bams_star; bam_merge as merge_bams_hisat2; bam_merge as merge_bams_minimap2 } from './modules/bam_merge.nf'
include { stringtie2 as stringtie2_short; stringtie2 as stringtie2_long } from './modules/stringtie2.nf'
include { strawberry } from './modules/strawberry.nf'
include { trinity_gg } from './modules/trinity_genome-guided.nf'
include { psiclass } from './modules/psiclass.nf'
include { tama } from './modules/tama.nf'
include { unagi } from './modules/unagi.nf'
include { portcullis } from './modules/portcullis.nf'
include { braker2 } from './modules/braker2.nf'
include { taco as taco_stringtie2_short; taco as taco_strawberry_short; taco as taco_stringtie_long } from './modules/taco.nf'
include { mikado_prepare } from './modules/mikado_prepare.nf'
include { prepare_mikado_file } from './modules/prepare_mikado_file.nf'



// Include QC tools
include { busco as busco_vir; busco as busco_emb; busco as busco_eud } from './modules/busco.nf'
include { quast } from './modules/quast.nf'
include { multiqc } from './modules/multiqc.nf'
include { kat } from './modules/kat.nf'



workflow {

/////////////// BASECALLING PIPELINE ///////////////
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
/////////////// CLEANUP PIPELINE ///////////////
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
/////////////// ASSEMBLY PIPELINE ///////////////
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

    wtdbg(params.fastq, params.genome_size)

    assembly = wtdbg.out.assembly
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

//////// LONG READ POLISHING ////////

  if ( params.polish ) {
    // racon and medaka polishing
    minimap2_1(params.fastq, assembly)
    racon1(params.fastq, minimap2_1.out.map, assembly)

    medaka(params.fastq, racon1.out.assembly)

    //medaka.out.assembly

    assembly = medaka.out.assembly

  } else if (params.racon_polish) {
    // only racon polishing
    minimap2_1(params.fastq, assembly)

    racon1(params.fastq, minimap2_1.out.map, assembly)

    assembly = racon1.out.assembly
  }else if (params.medaka_polish) {
    // only medaka polishing
    medaka(params.fastq, assembly)

    assembly = medaka.out.assembly
  }

//////// SHORT READ POLISHING ////////
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

        bam_coverage.out.coverage.subscribe { println "Short read coverage is ${it}x." }

        coverage = bam_coverage.out.coverage
      }
      else if ( params.short_polish_map == "minimap2" ) {
        // Mapping with minimap2 -ax sr

        minimap2_sr( short_r, assembly )

        short_bam = minimap2_sr.out.bam
        short_baidx = minimap2_sr.out.baidx

        bam_coverage( short_bam )

        bam_coverage.out.coverage.subscribe { println "Short read coverage is ${it}x." }

        coverage = bam_coverage.out.coverage

      }
      else if ( params.short_polish_map == "longranger" | params.short_polish_map == "10x" ) {
        // Longranger mapping for 10x data
        // Longranger is not distributed under conda or docker. Local installation?
        error 'Sorry. Longranger mapping is not yet implemented.'

      }
      else {
        error 'Unknown mapping method: ${params.short_polish_map}. Please choose from bwa, minimap2 or contact developers.'
      }

    }

    if ( params.short_polish == 'freebayes' | params.short_polish == 'vgp' | params.short_polish == true ) {
      // Verterae Genome Project polishing pipeline with freebayes and bcftools

      contig_index = Channel.from(1..100) // Split into 100 parallel processes

      freebayes_call( assembly, coverage.first(), short_bam.first(), short_baidx.first(), contig_index )

      // how collectFile works? Concatenate files or colelcts all files?
      //freebayes_call.out.collectFile(name: 'concat_list.txt', newLine: true, sort: true).set { bcf_list }

      freebayes_consensus( assembly, freebayes_call.out.collect() )

      polished_assembly = freebayes_consensus.out.assembly

    }
    else if ( params.short_polish == 'nextpolish' ) {
      // Nextpolish polishing
      if ( !params.nextpolish ) {
        error 'NextPolish source is not provided. Please install NextPolish locally and provide as --nextpolish /PATH/TO/nextpolish or consider using another polishing method (e.g. vgp).'
      }

      error 'Sorry. Nextpolish is not yet implemented.'

    }
    else if ( params.short_polish == 'pilon' ) {
      // pilon polishing

      pilon( assembly, short_bam, short_baidx )

      polished_assembly = pilon.out.assembly

    }
    else if ( params.short_polish == 'polca' ) {
      // POLCA polishing
      error 'Sorry. POLCA is not yet implemented.'
    }
    else if ( params.short_polish == 'hypo' ) {
      // HyPo polishing

      hypo( assembly, short_r, params.fastq, params.genome_size, short_bam, short_baidx, coverage )

      polished_assembly = hypo.out.assembly
    }

  }
  else {
    polished_assembly = assembly
  }


//////// ASSEMBLY QC ////////
  if ( !params.skip_qc) {
    // QC
    quast(polished_assembly)

    // busco_lineages = Channel.fromList = ['eudicots_odb10', 'embryophyta_odb10', 'viridiplantae_odb10']
    // busco(polished_assembly, busco_lineages, "genome")

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
/////////////// GENOME QC PIPELINE ///////////////
else if ( params.mode == 'genome_qc' ) {

  // Run quast and busco on an assembled genome
  quast(params.genome)
  busco_eud(params.genome, "eudicots_odb10", "genome")
  busco_emb(params.genome, "embryophyta_odb10", "genome")
  busco_vir(params.genome, "viridiplantae_odb10", "genome")

  multiqc(quast.out.summary.mix(busco_eud.out, busco_emb.out, busco_vir.out).collect(), "$baseDir/${params.outdir}")

  if ( params.short_reads ) {

    short_reads = Channel.fromFilePairs( params.short_reads )

    kat( params.genome, short_reads )
  }

}

/////////////// ANNOTATION PIPELINE ///////////////

else if ( params.mode == 'annotation' ) {
  log.info "Starting annotation protocol ... "

  // Genome file is provided
  if ( params.genome != false ) {
    // Genome file is provided, ru
    log.info "Genome file provided: ${params.genome}"
    log.info "Annotation pipeline is starting ..."

    // edta repeat masking
    if ( !params.skip_softmask ) {
      log.info "Soft masking repeats ... "
      edta(params.genome)

      masked_genome = edta.out.masked
      masked_gff = edta.out.te_anno
    }

  } else {
    log.info "No genome is provided."
    error "This is not yet supported."
  }

  if ( params.short_reads != false ) {
    // Map short reads to genome and assemble transcripts
    Channel
      .fromFilePairs( params.short_reads, size: params.single_end ? 1 : 2 )
      .ifEmpty { exit 1, "Reads are not provided correctly ${params.short_reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
      .set { short_reads }

    //fastp ?

    // STAR mapping
    star_idx(params.genome)
    star_align(star_idx.out, short_reads)

    //merge bam files to one bam file
    merge_bams_star( star_align.out.bam.collect(), "star" )
    // SHOULD NOT MERGE BAM, BUT ASSEMBLY INDEPENDENTLY AND USE TACO TO MERGE

    // HISAT2 mapping
    //hisat2_idx(params.genome)
    //hisat2_align(hisat2_idx.out, short_reads)

    //merge bam files to one bam file
    //merge_bams_hisat2( hisat2_align.out.bam.collect(), "hisat2" )

    // Transcript assemblies

    //stringtie2_short( params.genome, merge_bams_star.out.bam, '' )
    // ##### #TRY ##########
    stringtie2_short( params.genome, star_align.out.bam, '' ) // The last '' is a placeholder for long read settings
    taco_stringtie2_short( stringtie2_short.out.gtf.collect(), "stringtie2_short" )

    strawberry( params.genome, star_align.out.bam )
    taco_strawberry_short( strawberry.out.gtf.collect(), "strawberry_short" )

    trinity_gg( params.genome, merge_bams_star.out.bam )

    //psiclass( params.genome, star_align.out.bam.collect() )

    portcullis( params.genome, merge_bams_star.out.bam )


    // This should move out and merged with other gtf files from long reads and braker
    taco_stringtie2_short.out.gtf
      .collect()
      .mix( taco_strawberry_short.out.gtf.collect(), trinity_gg.out.gtf.collect() )
      .collect()
      .set { short_gtf }
    short_gtf.subscribe { println "Gene models generated from long reads:\n$it" }

    portcullis.out.junctions
      .set { junctions }

    // blast_mikado( "blastx", params.protein, mikado_prepare.out.fasta )

    // mikado_pick( params.genome, blast.mikado.out.xml, mikado_prepare.out.config, mikado_prepare.out.gtf, mikado_prepare.out.fasta )

    // PASA to generate gff3


    //outputs:
    // short.bam
    // short_transcripts.fasta
    // short_transcripts.gff3

  } else {
    Channel.from('')
      .set { short_gtf }
    Channel.from('')
      .set { junctions }
  }

  if ( params.ont_reads != false ) {

    //nanoq() ?
    //pychopper() ?
    Channel
      .fromFilePairs( params.ont_reads, size: 1 )
      .ifEmpty { exit 1, "ONT reads are not provided correctly ${params.ont_reads}\nNB: Path needs to be enclosed in quotes!" }
      .set { ont_reads }

    minimap_rna( params.genome, ont_reads )

    merge_bams_minimap2( minimap_rna.out.bam.collect(), "minimap2" )

    stringtie2_long( params.genome, minimap_rna.out.bam, '-L -l STRGL' )
    taco_stringtie_long( stringtie2_long.out.gtf.collect(), "stringtie2_long" )
    // tama( params.genome, minimap_rna.out.bam )
    unagi( params.genome, ont_reads )

    // Stringtie2 or Cufflinks or CLASS2 or UNAGI to reconstruct transcripts (fasta)
    // TACO or Mikado to merge transcript models

    // predict CDS? transdecoder()

    taco_stringtie_long.out.gtf
      .collect()
      .mix( unagi.out.gtf.collect() )
      .collect()
      .set { ont_gtf }
    ont_gtf.subscribe { println "Gene models generated from long reads:\n$it" }
  }  else {
      Channel.from('')
        .set { ont_gtf }
    }

  if ( !params.skip_abinitio ) {
    if ( ( params.short_reads || params.ont_reads ) && params.protein) {
      mark = "--etpmode"
    } else if ( params.protein ) {
      mark = "--epmode"
    } else if ( params.short_reads || params.ont_reads ) {
      mark = "--prg=gth --gth2traingenes"
    } else {
      mark = "--esmode --prg=gth --gth2traingenes"
    }
    // run BRAKER2
    if ( params.short_reads ) {
      braker2( params.genome, merge_bams_star.out.bam, mark )
    } else {
      braker2( params.genome, null, mark )
    }
    // braker2( params.genome )
    //outputs:
    // ab_initio.gff3

    braker2.out.gtf.set { braker_gtf }
    braker_gtf.subscribe { println "Gene models generated by BRAKER2:\n$it" }
  } else {
    Channel.from('')
      .set { braker_gtf }
  }


  braker_gtf
    .collect()
    .mix(short_gtf, ont_gtf)
    .collect()
    .set { all_gtf }

  prepare_mikado_file( all_gtf )

  mikado_prepare( params.genome, all_gtf, prepare_mikado_file.out.gtf_list, params.mikado_scoring, junctions )

  // transdecoder + blastx
  //mikado_pick()

  // Update annotation with funannotation or PASA for UTRs and isoforms
  // Use merged transcripts + cDNA evidence from long reads

  // Predict CDS?

  // Functional annotation?
  // interproscan()
  // eggnog()

  // trnascan-se()

  //small RNA
  // https://github.com/MikeAxtell/ShortStack
  // https://github.com/miRkwood-RNA/miRkwood
  if ( params.srna ) {
    if ( params.srna_reads ) {
      // Annotation based on sequencing data

    } else {
      // Ab initio annotation of sRNAs

    }
  }


}

/////////////// Errors when no or invalid mode is provided. ///////////////

else if ( !params.mode ) {
  error 'Running mode is not supplied. Please specify with --mode'
} else {
  error 'Invalid running method: $params.mode. Currently supported modes: basecalling, assembly, annotation, expression.'
}


}


workflow.onComplete {
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> Nanopore $params.mode pipeline finished SUCCESSFULLY after $workflow.duration ."
      log.info "[$workflow.complete] >> "
      log.info "[$workflow.complete] >> You can find further help on https://github.com/palfalvi/nanoporeseq"
    } else {
      log.info "[$workflow.complete] >> The script quit with ERROR after ${workflow.duration}."
      log.info "[$workflow.complete] >> Please revise your code and resubmit jobs with the -resume option or reach out for help at https://github.com/palfalvi/nanoporeseq."
    }
}
