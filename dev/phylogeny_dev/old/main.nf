#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// include functions
include {printHelp} from "${params.nfpath}/modules/help.nf"

// import subworkflows
include {bacteria_phylogeny} from "${params.nfpath}/workflows/workflow.nf"

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
    //Step1- create a Channel based on content of replicons.tsv
    replicons_files = Channel.fromPath(params.result+"/replicons.tsv", checkIfExists:true).splitCsv(sep:'\t', header: true)
    replicons_ch = replicons_files.map { row -> tuple(row.replicon, row.Sample, row.fasta_file) }
    //replicons_ch.groupTuple(by: 0).set{ replicons_ch }
    //replicons_ch.view()
    
    main:
    //Step2- launch the appropriate workflow
    if ( params.workflow == 'bacteria_phylogeny') {
        bacteria_phylogeny(replicons_ch)
    } 
    
   
}