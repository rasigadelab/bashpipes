def printHelp() {

  log.info"""
  Usage:
    nextflow run main.nf -profile (standard) [workflow-options]

  Description:
    Nextflow pipeline for de novo assembly of Nanopore/Illumina reads. 
    Specific to microbial isolates.
    Included one workflow, producing Nanopore/Illumina hybrid assembly,
    with Flye assembling Nanopore reads and Pilon correcting with Illumina reads.

  Nextflow arguments:
    

  Mandatory workflow arguments:


  Variant workflow options:
    Mandatory:

    Optional:

  """.stripIndent()
}