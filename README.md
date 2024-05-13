# Bioinformatic pipelines for epidemiological study of clinical bacterial isolates

Common pipelines currently used and developed in PHE3ID team. The team focus on prokarytic studies in a clinical context of antibiotic resistance emergence and surveillance.

Two types of analyses can be led with these pipelines.
* Epitrack pipeline: _de novo_ assembly of prokaryotic genome using Illumina short-read data, followed by global genome annotations (taxonomy ; ARGs identification ; plasmid identification)
* Resistrack pipeline: hybrid _de novo_ assembly of prokaryotic genome using Nanopore long-read data and corrected by Illumina data, followed by global genome annotations (taxonomy ; ARGs identification ; plasmid identification)

Some features under development...
*  Phylogenetic analyses: phylogenetic tree based on Maximum-Likelihood reconstruction (IQtree) and core SNP analyses (Snippy).
*  Automated reports: reports to properly present results of phylogenetic analyses, applied to clonal investigation among a batch of clinical samples.

All pipelines are running with Nextflow and Singularity containers. More details on specifications here after.

## Epitrack pipeline

![diagramme_epitrack](https://github.com/rasigadelab/bashpipes/assets/120658937/d8d638e7-da21-4144-b3ce-1c248822fd40)


## Resistrack pipeline

## Phylogeny pipeline
