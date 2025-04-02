#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// include functions
include {printHelp} from "${params.nfpath}/modules/help.nf"
include {make_sample_dir} from "${params.nfpath}/modules/util.nf"

// import subworkflows
include {plasmid_compa} from "${params.nfpath}/workflows/workflow.nf"

def raiseError ( value ) {
    sleep(2000)
    println(value)
    System.exit(1)
}

if (params.help) {
    printHelp()
    exit 0
}

// main workflow
workflow {

    main:

    //Step1- create a Channel based on content of plasmid_locations.tsv
    plasmid_locations = Channel.fromPath(params.result+"/plasmid_locations.tsv", checkIfExists:true).splitCsv(sep:'\t', header: false)
    //plasmid_locations.groupTuple(by: [0,1,3]).set{ plasmid_ch }
    plasmid_compa(plasmid_locations)
     
   
}