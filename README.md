# Bioinformatic pipelines for epidemiological study of clinical bacterial isolates

Pipelines currently used and developed in PHE3ID team. The team focus on prokaryotic studies in a clinical context of antibiotic resistance emergence and surveillance.

Two types of analyses can be led with these pipelines.
* Illumina-only pipeline: _de novo_ assembly of prokaryotic genome using Illumina short-read data, followed by global genome annotations (taxonomy ; ARGs identification ; plasmid identification)
* Nano-Illumina pipeline: hybrid _de novo_ assembly of prokaryotic genome using Nanopore long-read data and corrected by Illumina data, followed by global genome annotations (taxonomy ; ARGs identification ; plasmid identification)

Some features under development...
*  Phylogenetic analyses: phylogenetic tree based on Maximum-Likelihood reconstruction (IQtree) and core SNP analyses (Snippy).
*  Automated reports: reports to properly present results of phylogenetic analyses, applied to clonal investigation among a batch of clinical samples.
*  Resistrack pipeline: Adding a module for Nanopore long read trimming with Porechop.

All pipelines are running with Nextflow and Singularity containers. More details on specifications here after.

## Illumina-only pipeline

Illumina-only pipeline has been developed to analyze Illumina short-read data and to assemble reads in order to get an assembly of the studied genome. 
First, quality of short reads are assessed using fastp (https://github.com/OpenGene/fastp, version 0.23.4) and adapters trimmed with Trimmomatic (https://github.com/usadellab/Trimmomatic, version 0.39). Reads are then assembled with SPAdes assembler, with option « isolate » (https://github.com/ablab/spades, version 3.15.5). Contigs with less than 500bp are removed from analysis. Afterwards, contigs are circularised with Circlator (https://github.com/sanger-pathogens/circlator, version 1.5.5). Quality metrics (N50, number of contigs,…) on the final assembly are computed with Quast (https://github.com/ablab/quast, version 5.2.0). 

![diagramme_epitrack](https://github.com/rasigadelab/bashpipes/assets/120658937/53c2ba66-5a22-4447-b2d3-d7abd37cb177)

Genomes are annonated to retrieve the taxonomy of the sample using Sourmash (https://github.com/sourmash-bio/sourmash, version 4.8.5). Genomes are typed with classical Multi-Locus Sequence Typing (MLST - https://github.com/tseemann/mlst, version 2.23.0). Antibiotic resistance genes (ARGs) are detected with NCBI AMRfinder+ (https://github.com/ncbi/amr, version 3.12.8). Global genome annotations are computed with Bakta (https://github.com/oschwengers/bakta, version 1.9.2). Finally, contigs corresponding to plasmids are identified using Mob-Recon of the Mob-Suite (https://github.com/phac-nml/mob-suite, version 3.1.8). 

## Nano-Illumina pipeline

Nano-Illumina pipeline has been developed to build genome assemblies based on Nanopore and Illumina sequencing data.
Quality of short reads are assessed using fastp (https://github.com/OpenGene/fastp, version 0.23.4) and adapters trimmed with Trimmomatic (https://github.com/usadellab/Trimmomatic, version 0.39). Long reads are filtered with Filtlong (https://github.com/rrwick/Filtlong, version 0.2.1) to keep only the 80% bases with best quality and reads longer than 5kb. Statistics on filtered reads are computed with NanoPlot (https://github.com/wdecoster/NanoPlot, version 1.42.0). 
Nanopore long reads are then assembled with Flye assembler (https://github.com/fenderglass/Flye, version 2.9.1) to get a first draft assembly. « --meta » option was used to accomodate the assembler to variation of coverage between reads of plasmids (shorter and thus potentially more present in the dataset) and reads of chromosome (longer and thus potentially less present in the dataset). Short Illumina reads are mapped back on Nanopore assembly with Bowtie2 (https://github.com/BenLangmead/bowtie2, version 2.2.5) and the draft assembly is corrected based on short reads with Pilon software (https://github.com/broadinstitute/pilon, version 1.24). Afterwards, contigs are circularised with Circlator (https://github.com/sanger-pathogens/circlator, version 1.5.5).

![nano_pipe](https://github.com/rasigadelab/bashpipes/assets/120658937/c636cc96-b2ab-483c-9461-21762b27ca74)

The same global genomic annotations are made as for Epitrack pipeline. Taxonomy is retrieved with Sourmash. Each genome is typed by MLST. ARGs are detected with AMRFinder+ and genes are annotated with Bakta. Contigs corresponding to plasmids and chromosome are identified with Mob-Recon.

## Phylogeny pipeline [under development]

Samples that are potentially clonal are analyzed together to compare the number of SNPs localized in their chromosomes. Samples are first clustered in small groups of genetically near samples. Groups are created by computing the mash distance (https://github.com/marbl/Mash, version 2.3) between each isolates pair in the cluster and by applying a single-linkage clustering with a cut at a specific threshold (in most cases 0.005). This enables to identify outliers or samples with an assembly of bad quality. For each created subcluster a reference sample is chosen : it is defined as the isolate the nearest from the most of the other samples in the subcluster. 

Reads of reference samples are mapped back to reference genome assemblies with Bowtie2 (https://github.com/BenLangmead/bowtie2, version 2.2.5) and genome assemblies are corrected based on alignments with Pilon (https://github.com/broadinstitute/pilon, version 1.24). Repeat regions are identified on reference sample assembly with RepeatMasker (https://github.com/rmhubley/RepeatMasker, version 4.1.5). Then, SNP calling and core SNPs detection are made on the batch of samples of each subcluster, taking the chosen reference sample as reference genome for SNP calling. This step is computed with Snippy (https://github.com/tseemann/snippy, version 4.6.0) and results in a full SNP alignment of the batch of samples. Clonal samples are defined as samples distant from less than 15 SNPs.

Samples of the main cluster are then analyzed together. SNP calling is made with Snippy on the whole batch of genomes, with the reference sample of the biggest minicluster used as reference genome. SNPs alignment is then processed by Iqtree (https://github.com/iqtree/iqtree2, version 2.2.0.3) to reconstruct phylogeny of the batch of samples. GTR model is used as parameter for tree search. Phylogenetic tree is corrected from recombination with ClonalFrameML (https://github.com/xavierdidelot/clonalframeml, version 1.12). A datation estimation is computed at the end of the pipeline by treetime (https://github.com/neherlab/treetime, version 0.11.1).

![phylo_pipe](https://github.com/rasigadelab/bashpipes/assets/120658937/7c464a6a-13f0-4945-9aba-44f9a28989f9)

