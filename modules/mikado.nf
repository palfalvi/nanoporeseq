include { mikado_prepare } from '../modules/mikado_prepare.nf'
include { mikado_serialise } from '../modules/mikado_serialise.nf'
include { blast_makedb } from '../modules/blast_makedb.nf'
include { blast } from '../modules/blast.nf'
include { transdecoder } from '../modules/transdecoder.nf'
include { mikado_pick } from '../modules/mikado_pick.nf'

workflow mikado {
    take:
      genome
      gffs
      scoring
      //junction

    main:
      mikado_prepare(  genome, gffs, scoring  )

      if ( params.protein ) {

        blast_makedb( params.protein, 'prot' )
        blast( blast_makedb.out, mikado_prepare.out.fasta.splitFasta( by: 5000, file: true ), 'blastx' )
        blast.out.blast.collectFile(name: 'mikado_prepared.blast.tsv', newLine: true).set { blastp }

      } else {

        Channel.from([]).set { blastp }

      }

      transdecoder( mikado_prepare.out.fasta )

      mikado_serialise( genome,  blastp, transdecoder.out.bed, scoring, mikado_prepare.out.mikado )
      mikado_pick( genome, mikado_serialise.out.mikado, mikado_prepare.out.gtf, scoring,   )

    emit:
        loci = mikado_pick.out.loci
        subloci = mikado_pick.out.subloci
        metrics = mikado_pick.out.metrics
        scores = mikado_pick.out.scores
}
