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

        mikado_blastp_makeref( params.protein, 'prot' )
        mikado_blastp( mikado_prepare.out.fasta.splitFasta( by: 5000, file: true ), mikado_blastp_makeref.out, 'blastx' )
        mikado_blastp.out.blast.collectFile(name: 'mikado_prepared.blast.tsv', newLine: true).set { blastp }

      } else {

        Channel.from([]).set { blastp }

      }

      mikado_transdecoder( mikado_prepare.out.fasta )

      mikado_serialise( genome, mikado_prepare.out.config, blastp, mikado_transdecoder.out.bed, scoring )
      mikado_pick( genome, mikado_prepare.out.config, mikado_serialise.out.db, mikado_prepare.out.gtf, scoring  )

    emit:
        loci = mikado_pick.out.loci
        subloci = mikado_pick.out.subloci
        metrics = mikado_pick.out.metrics
        scores = mikado_pick.out.scores
}
