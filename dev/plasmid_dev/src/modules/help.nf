/*
 * Help utilities for Bacteria Plasmid Comparison Pipeline
 *
 * Description:
 *   Defines functions to display usage and CLI help messages.
 *
 * License:
 *   AGPL-3.0-only
 */

def printHelp() {

  log.info"""
  ${workflow.manifest.name ?: 'Pipeline'}
  
  Usage
    nextflow -C nextflow.config run main.nf -params-file params_plasmid_compa.json -profile <profile> [options]

  Description:
    ${workflow.manifest.description ?: 'Plasmid comparison pipeline based on prokaryotic reads'}

  Required arguments:
    -params-file        Path to parameters configuration
    -C                  Path to nextflow configuration file

  Optional arguments:
    --help         Display this help message
    --nfpath       Path to pipeline resources

  Profiles:
    standard       Default execution
    developer      Run with low resources consumption (CPU-RAM)

  """.stripIndent()
}
