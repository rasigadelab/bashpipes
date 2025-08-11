#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {quality_fastp} from "${params.nfpath}/modules/module.nf"
include {trim_trimmomatic} from "${params.nfpath}/modules/module.nf"
include {resync_bbmap} from "${params.nfpath}/modules/module.nf"
include {assembly_spades} from "${params.nfpath}/modules/module.nf"
include {filter_contigs_bbmap} from "${params.nfpath}/modules/module.nf"
include {qc_quast} from "${params.nfpath}/modules/module.nf"
include {fixstart_circlator} from "${params.nfpath}/modules/module.nf"

include {mlst_sequence_typing} from "${params.nfpath}/modules/module.nf"
include {classify_sourmash} from "${params.nfpath}/modules/module.nf"
include {amr_typer_amrfinder} from "${params.nfpath}/modules/module.nf"
include {annotate_bakta} from "${params.nfpath}/modules/module.nf"
include {mge_mob_recon} from "${params.nfpath}/modules/module.nf"

// workflow script
workflow bacteria_denovo {
    take:
        ch_illumina
        ch_ont

    main:
       
       //Step0- Reads quality and fastq trimming
       if ( params.fastp ) {
            quality_fastp(ch_illumina)
            quality_fastp.out.illumina_reads.set{ ch_illumina }
       }

       if ( params.trimming ) {
            trim_trimmomatic(ch_illumina)
            trim_trimmomatic.out.illumina_trimmed.set{ ch_illumina }
            resync_bbmap(ch_illumina)
            resync_bbmap.out.illumina_resync.set{ ch_illumina }
       }

        //StepA- De novo assembly
        if ( params.assembler == 'spades' ) {
            assembly_spades(ch_illumina)
            assembly_spades.out.draft_assembly.set{ ch_denovo_assembly }
            filter_contigs_bbmap(ch_denovo_assembly)
            filter_contigs_bbmap.out.draft_assembly.set{ ch_denovo_assembly }
        }

        //StepB- Assembly QC
        if ( params.quast ) {
            qc_quast(ch_denovo_assembly)
            qc_quast.out.draft_assembly.set{ ch_denovo_assembly }
        }
        
        //StepC- Circularization of contigs
        if ( params.circularization ) {
            fixstart_circlator(ch_denovo_assembly)
            fixstart_circlator.out.realigned_assembly.set{ ch_final_assembly }
        }
        
        //StepD- MLST Sequence Typing
        if ( params.mlst ) {
            mlst_sequence_typing(ch_final_assembly)
            mlst_sequence_typing.out.final_assembly.set{ ch_final_assembly }
        }
        
        //StepE- Taxon classification
        if ( params.sourmash ) {
            classify_sourmash(ch_final_assembly)
            classify_sourmash.out.sample_taxonomy.set{ ch_sample_taxonomy }
        }

        //StepF- AMR genes annotation (+other genes of interest)
        if ( params.amr_typer ) {
            amr_typer_amrfinder(ch_final_assembly.join(ch_sample_taxonomy))
            amr_typer_amrfinder.out.final_assembly.set{ ch_final_assembly }
        }
        
        //StepG- Genome annotation
        if ( params.annotate ) {
            annotate_bakta(ch_final_assembly)
            annotate_bakta.out.final_assembly.set{ ch_final_assembly }
        }

        //StepH- MGE Analysis (mobile genetic elements)
        if ( params.mge ) {
            mge_mob_recon(ch_final_assembly)
            mge_mob_recon.out.samples_list.set{ ch_samples }
        }
}