def printHelp() {

  log.info"""
  Usage:
    nextflow run main.nf -profile (standard) [workflow-options]

  Description:
    Nextflow pipeline for phylogeny of Nanopore/Illumina reads. 
    Specific to prokaryotes.

  Nextflow arguments:
    

  Mandatory workflow arguments:


  Variant workflow options:
    Mandatory:

    Optional:

  """.stripIndent()
}