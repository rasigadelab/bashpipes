#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {assembler_flye} from "${params.nfpath}/modules/module.nf"

// workflow script
workflow bacteria_denovo {
    take:
        ch_illumina
        ch_ont

    main:


        ch_illumina.map{ it[0] }.set{ ch_sampleId }
    /*
        if ( params.assembler == 'flye' ) {
            assembler_flye(ch_ont, ch_sampleId)
            assembler_flye.out.set{ ch_ont }
        } else if ( params.assembler == 'spades' ) {

        }
    */
}