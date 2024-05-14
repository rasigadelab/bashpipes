def printHelp() {

  log.info"""
  Usage:
    nextflow run main.nf -profile (standard) [workflow-options]

  Description:
    Nextflow pipeline for de novo assembly of Illumina reads. 
    Specific for bacterial isolates.
    Included one workflow, producing Illumina-only assembly with SPAdes assembler.

  Nextflow arguments:
    

  Mandatory workflow arguments:


  Variant workflow options:
    Mandatory:

    Optional:

  """.stripIndent()
}
