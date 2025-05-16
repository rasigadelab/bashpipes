#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules
include {gather_plasmid_seq} from "${params.nfpath}/modules/module.nf"
include {align_to_reference} from "${params.nfpath}/modules/module.nf"
include {visualize_circos} from "${params.nfpath}/modules/module.nf"
include {change_circos_config} from "${params.nfpath}/modules/module.nf"

// workflow script
workflow plasmid_compa {
    take:
        ch_plasmids

    main:

       // Step0- Gathering of plasmid sequences
       if ( params.gathering ) {
            gather_plasmid_seq(ch_plasmids)
            gather_plasmid_seq.out.plasmid_seq.set{ ch_plasmids }
            ch_plasmids.groupTuple(by: [0,1])
                       .flatMap { key1, key2, values -> values.collate(8).indexed().collect { idx, chunk -> [key1, key2, chunk, "Batch${idx + 1}"]}}
                       .set{ ch_plasmids }
       }

       // Step1- Alignement à la référence
        if ( params.mapping ) {
            align_to_reference(ch_plasmids)
            align_to_reference.out.plasmid_seq.set { ch_plasmids }
        }

        // Step2- Circos visualisation
        if ( params.visual ) {
            ch_plasmids.flatMap { key1, key2, values -> values.flatten().collect { chunk -> [key1, key2, chunk]}}
                       .groupTuple(by: [0,1])
                       .set{ ch_plasmids }
            visualize_circos(ch_plasmids)
            visualize_circos.out.plasmids_info.set{ch_plasmids}
            change_circos_config(ch_plasmids)
        }

            

}       