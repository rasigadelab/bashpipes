#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {assembly_flye} from "${params.nfpath}/modules/module.nf"
include {assembly_spades} from "${params.nfpath}/modules/module.nf"
include {map_bowtie2} from "${params.nfpath}/modules/module.nf"
include {polish_pilon} from "${params.nfpath}/modules/module.nf"
include {qc_quast} from "${params.nfpath}/modules/module.nf"
include {fixstart_circlator} from "${params.nfpath}/modules/module.nf"


// workflow script
workflow bacteria_denovo {
    take:
        ch_illumina
        ch_ont

    main:
       
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

    

        
    
}