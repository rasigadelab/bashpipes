#!/home/afischer/tools/nextflow

// enable dsl2
nextflow.enable.dsl = 2

// include functions
include {printHelp} from "${params.nfpath}/modules/help.nf"
include {make_sample_dir} from "${params.nfpath}/modules/util.nf"

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
    /*
    //Step2- launch process to create a directory for each sample in /genomes folder
    make_sample_dir(fastq_ch)
    //Step3- create a Channel [sample, ONT/R1/R2_output_file_path]
    make_sample_dir.out.set{ new_ch }
    //Step3- create a Channel for each type of reads (ONT, R1 or R2)
    illumina_R1 = new_ch.unique().filter{ it[1] =~/.*\_R1.fastq.gz$/ }
    illumina_R2 = new_ch.unique().filter{ it[1] =~/.*\_R2.fastq.gz$/ }
    ont_ch = new_ch.unique().filter{ it[1] =~/.*\_ONT.fastq.gz$/ }
    //Step4- create a final Channel joining R1 and R2 Illumina reads together by sample.
    illumina_ch = illumina_R1.join(illumina_R2)
    */
    
    main:
    //Step5- launch the appropriate workflow
    if ( params.workflow == 'bacteria_phylogeny') {
        bacteria_phylogeny(replicons_ch)
    } 
    
   
}