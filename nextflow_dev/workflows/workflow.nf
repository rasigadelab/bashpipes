#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {assembly_flye} from "${params.nfpath}/modules/module.nf"
include {assembly_spades} from "${params.nfpath}/modules/module.nf"


// workflow script
workflow bacteria_denovo {
    take:
        ch_illumina
        ch_ont

    main:
        //ch_illumina.map{ it[0] }.set{ ch_sampleId }
       
    
        if ( params.assembler == 'flye' ) {
            assembly_flye(ch_ont)
            assembly_flye.out.draft_assembly.set{ ch_draft_assembly }
        } else if ( params.assembler == 'spades' ) {
            assembly_spades(ch_illumina)
            assembly_spades.out.draft_assembly.set{ ch_draft_assembly }
        }
    
}