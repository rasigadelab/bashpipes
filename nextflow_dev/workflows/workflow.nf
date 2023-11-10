#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {filter_filtlong} from "${params.nfpath}/modules/module.nf"
include {assembly_flye} from "${params.nfpath}/modules/module.nf"
include {assembly_spades} from "${params.nfpath}/modules/module.nf"
include {map_bowtie2} from "${params.nfpath}/modules/module.nf"
include {polish_pilon} from "${params.nfpath}/modules/module.nf"
include {qc_quast} from "${params.nfpath}/modules/module.nf"
include {fixstart_circlator} from "${params.nfpath}/modules/module.nf"

include {mlst_sequence_typing} from "${params.nfpath}/modules/module.nf"
include {classify_sourmash} from "${params.nfpath}/modules/module.nf"
include {amr_typer_amrfinder} from "${params.nfpath}/modules/module.nf"
include {annotate_prokka} from "${params.nfpath}/modules/module.nf"
include {mge_mob_recon} from "${params.nfpath}/modules/module.nf"
include {quality_fastqc} from "${params.nfpath}/modules/module.nf"
include {trim_trimmomatic} from "${params.nfpath}/modules/module.nf"
include {stats_nanoplot} from "${params.nfpath}/modules/module.nf"

include{cleaning} from "${params.nfpath}/modules/module.nf"

// workflow script
workflow bacteria_denovo {
    take:
        ch_illumina
        ch_ont

    main:
       
       //Step0- Reads quality and fastq trimming
       if ( params.fastqc ) {
            quality_fastqc(ch_illumina)
            quality_fastqc.out.illumina_reads.set{ ch_illumina }
       }

       if ( params.trimming ) {
            trim_trimmomatic(ch_illumina)
            trim_trimmomatic.out.illumina_trimmed.set{ ch_illumina }
       }

       if ( params.nano_filtering ) {
            filter_filtlong(ch_ont)
            filter_filtlong.out.filtered_nanopore.set{ ch_ont }
       }

       if ( params.nanoplot ) {
            stats_nanoplot(ch_ont)
            stats_nanoplot.out.nanopore_reads.set{ ch_ont }
       }

        //StepA- De novo assembly
        if ( params.assembler == 'flye' ) {
            assembly_flye(ch_ont)
            assembly_flye.out.draft_assembly.set{ ch_draft_assembly }
        } else if ( params.assembler == 'spades' ) {
            assembly_spades(ch_illumina)
            assembly_spades.out.draft_assembly.set{ ch_denovo_assembly }
        }

        //Case: Nanopore-Illumina hybrid assembly
        if ( params.hybrid_assembly ) {
            //StepB- Illumina mapping on ONT draft assembly
            if ( params.map ) {
                map_bowtie2( ch_draft_assembly.join(ch_illumina) )
                map_bowtie2.out.sorted_bam_files.set{ ch_sorted_bam_files }
            }

            //StepC- Polishing draft assembly with Illumina reads
            if ( params.polish ) {
                polish_pilon(ch_draft_assembly.join(ch_sorted_bam_files))
                polish_pilon.out.polished_assembly.set{ ch_denovo_assembly }
            }
        //Case: Illumina-only assembly
        } else {
            //StepB- Assembly QC
            if ( params.quast ) {
                qc_quast(ch_denovo_assembly)
            }
        }
        
        //StepD- Circularization of contigs
        if ( params.circularization ) {
            fixstart_circlator(ch_denovo_assembly)
            fixstart_circlator.out.realigned_assembly.set{ ch_final_assembly }
        }
        
        //StepE- MLST Sequence Typing
        if ( params.mlst ) {
            mlst_sequence_typing(ch_final_assembly)
            mlst_sequence_typing.out.final_assembly.set{ ch_final_assembly }
        }
        
        //StepF- Taxon classification
        if ( params.sourmash ) {
            classify_sourmash(ch_final_assembly)
            classify_sourmash.out.sample_taxonomy.set{ ch_sample_taxonomy }
        }

        //StepG- AMR genes annotation (+other genes of interest)
        if ( params.amr_typer ) {
            amr_typer_amrfinder(ch_final_assembly.join(ch_sample_taxonomy))
            amr_typer_amrfinder.out.final_assembly.set{ ch_final_assembly }
        }
        
        //StepH- Genome annotation
        if ( params.annotate ) {
            annotate_prokka(ch_final_assembly)
            annotate_prokka.out.final_assembly.set{ ch_final_assembly }
        }

        //StepI- MGE Analysis (mobile genetic elements)
        if ( params.mge ) {
            mge_mob_recon(ch_final_assembly)
            mge_mob_recon.out.samples_list.set{ ch_samples }
        }

        //StepJ- Cleaning
        if ( params.clean ) {
            cleaning(ch_samples)
        }
    
}