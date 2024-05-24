def printHelp() {

  log.info"""
  Usage:
    nextflow run main.nf -profile (standard) [workflow-options]

  Description:
    Nextflow pipeline for de novo assembly of Nanopore/Illumina reads. 
    Specific to prokaryotes.
    Included two workflows, one producing Nanopore/Illumina hybrid assembly.
    Another one producing Illumina-only assembly.

  Nextflow arguments:
    

  Mandatory workflow arguments:


  Variant workflow options:
    Mandatory:

    Optional:

  """.stripIndent()
}