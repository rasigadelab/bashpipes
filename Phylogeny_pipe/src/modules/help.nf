/*
 * Help utilities for Bacteria Phylogeny Pipeline
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
  
  Usage A) MASH clustering
    nextflow -C nextflow.config run main.nf -params-file params_mash.json -profile <profile> [options]
  Usage B) Variant calling 
    nextflow -C nextflow.config run main.nf -params-file params_variant_calling.json -profile <profile> [options]
  Usage C) Phylogeny pipeline
    nextflow -C nextflow.config run main.nf -params-file params_phylogeny.json -profile <profile> [options]

  Description:
    ${workflow.manifest.description ?: 'Phylogeny pipeline for prokaryotic reads'}

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