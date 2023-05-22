#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {annotate_prokka} from "${params.nfpath}/modules/module.nf"


// workflow script
workflow bacteria_phylogeny {
    take:
        ch_replicons

    main:
       
       //Step0- Assembly annotation
       if ( params.annotate ) {
            annotate_prokka(ch_replicons)
            //annotate_prokka.out.illumina_reads.set{ ch_replicons }
       }
    
}