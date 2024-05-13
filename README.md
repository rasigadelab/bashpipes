# Bioinformatic pipelines for epidemiological study of clinical bacterial isolates

Common pipelines currently used and developed in PHE3ID team. The team focus on prokaryotic studies in a clinical context of antibiotic resistance emergence and surveillance.

Two types of analyses can be led with these pipelines.
* Epitrack pipeline: _de novo_ assembly of prokaryotic genome using Illumina short-read data, followed by global genome annotations (taxonomy ; ARGs identification ; plasmid identification)
* Resistrack pipeline: hybrid _de novo_ assembly of prokaryotic genome using Nanopore long-read data and corrected by Illumina data, followed by global genome annotations (taxonomy ; ARGs identification ; plasmid identification)

Some features under development...
*  Phylogenetic analyses: phylogenetic tree based on Maximum-Likelihood reconstruction (IQtree) and core SNP analyses (Snippy).
*  Automated reports: reports to properly present results of phylogenetic analyses, applied to clonal investigation among a batch of clinical samples.

All pipelines are running with Nextflow and Singularity containers. More details on specifications here after.

## Epitrack pipeline

Epitrack pipeline has been developed to analyze Illumina short-read data and to assemble reads in order to get an assembly of the studied genome. 
First, quality of short reads are assessed using fastp (https://github.com/OpenGene/fastp, version 0.23.4) and adapters trimmed with Trimmomatic (https://github.com/usadellab/Trimmomatic, version 0.39). Reads are then assembled with SPAdes assembler, with option « isolate » (https://github.com/ablab/spades, version 3.15.5). Contigs with less than 500bp are removed from analysis. Afterwards, contigs are circularised with Circlator (https://github.com/sanger-pathogens/circlator, version 1.5.5). Quality metrics (N50, number of contigs,…) on the final assembly are computed with Quast (https://github.com/ablab/quast, version 5.2.0). 

![diagramme_epitrack](https://github.com/rasigadelab/bashpipes/assets/120658937/53c2ba66-5a22-4447-b2d3-d7abd37cb177)

Genomes are annonated to retrieve the taxonomy of the sample using Sourmash (https://github.com/sourmash-bio/sourmash, version 4.8.5). Genomes are typed with classical Multi-Locus Sequence Typing (MLST - https://github.com/tseemann/mlst, version 2.23.0). Anti-microbial genes (ARGs) are detected with NCBI AMRfinder+ (https://github.com/ncbi/amr, version 3.12.8). Global genome annotations are computed with Bakta (https://github.com/oschwengers/bakta, version 1.9.2). Finally, contigs corresponding to plasmids are identified using Mob-Recon of the Mob-Suite (https://github.com/phac-nml/mob-suite, version 3.1.8). 

## Resistrack pipeline

## Phylogeny pipeline
