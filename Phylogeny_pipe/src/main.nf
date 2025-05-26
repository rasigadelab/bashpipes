#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// include functions
include {printHelp} from "${params.nfpath}/modules/help.nf"

// import subworkflows
include {bacteria_mash_clustering} from "${params.nfpath}/workflows/workflow.nf"
include {bacteria_variant_calling} from "${params.nfpath}/workflows/workflow.nf"
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

    main:
    
    if ( params.workflow == 'bacteria_mash_clustering') {
        //Step1- create a Channel based on content of replicons.tsv
        replicons_files = Channel.fromPath(params.result+"/replicons.tsv", checkIfExists:true).splitCsv(sep:'\t', header: true)
        replicons_ch = replicons_files.map { row -> tuple(row.replicon, row.Sample, row.fasta_file) }
        //Step2- launch the appropriate workflow
        bacteria_mash_clustering(replicons_ch)
    } else if ( params.workflow == 'bacteria_variant_calling') {
        //Step1- create a Channel based on content of replicons.tsv
        replicons_files = Channel.fromPath(params.result+"/replicons.tsv", checkIfExists:true).splitCsv(sep:'\t', header: true)
        replicons_ch = replicons_files.map { row -> tuple(row.Sample, row.replicon) }
        //Step2- Create a Channel based on content of mini_clusters.tsv
        clusters_files = Channel.fromPath(params.result+"/phylogeny/*/mash/mini_clusters.tsv", checkIfExists:true).splitCsv(sep:'\t', header: true)
        clusters_ch = clusters_files.map { row -> tuple(row.sample_names, row.mini_clusters, row.reference) }
        //Step3- Join both Channels
        input_ch = replicons_ch.join(clusters_ch).groupTuple(by: [1,2,3])
        //Step4- Launch appropriate workflow
        bacteria_variant_calling(input_ch)
    } else if ( params.workflow == 'bacteria_phylogeny') {
        //Step1- create a Channel based on content of replicons.tsv
        replicons_files = Channel.fromPath(params.result+"/replicons.tsv", checkIfExists:true).splitCsv(sep:'\t', header: true)
        replicons_ch = replicons_files.map { row -> tuple(row.replicon, row.Sample) }
        //Step2- Create a Channel based on content of mini_clusters.tsv
        reference_files = Channel.fromPath(params.result+"/phylogeny/*/mash/reference_for_phylogeny.tsv", checkIfExists:true).splitCsv(sep:'\t', header: true)
        reference_ch = reference_files.map { row -> tuple(row.replicon, row.reference) }
        //Step3- Gather samples for each replicon
        replicons_ch = replicons_ch.combine(reference_ch, by: 0).groupTuple(by: [0,2])
        //Step3- launch the appropriate workflow
        bacteria_phylogeny(replicons_ch)
    
    }
    
   
}