#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// include functions
include {printHelp} from "${params.nfpath}/modules/help.nf"
include {make_sample_dir} from "${params.nfpath}/modules/util.nf"

// import subworkflows
include {bacteria_denovo} from "${params.nfpath}/workflows/workflow.nf"

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
    fastq_locations = Channel.fromPath(params.result+"/Files_location.tsv", checkIfExists:true).splitCsv(sep:'\t', header: true)
    make_sample_dir(fastq_locations, params.result)
    
    if ( params.hybrid_assembly == true) {
        ont_path = "$params.result/genomes/*/*_ONT.fastq.gz"
        ont_ch = Channel.fromPath(ont_path, checkIfExists: true)
    } 
    else {
        ont_ch = Channel.of("0")
    }
    illumina_path = "$params.result/genomes/*/*_{R1,R2}.fastq.gz"
    illumina_ch = Channel.fromFilePairs(illumina_path, checkIfExists: true)
    
    main:
    if ( params.workflow == 'bacteria_denovo') {
        bacteria_denovo(illumina_ch, ont_ch)
    } 
    
}