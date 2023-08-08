#!/bin/bash
# End of PHYLOGENY Pipeline : After running Prokka on every samples
# Steps implemented in this script :
# 0- Pan-genome analysis
# 1- FASTA of Core genome consensus
# 2- Snippy input
# 3- Running Snippy
# 4- Running IQtree

# To use this script,
# Make sure you have 'core_genome.py' and 'snippy_multi_list.py' in same folder as this script
# Then type : ./end_phylogeny <replicon_name> <work_dir> [for example, ./end_phylogeny paeruginosa /mnt/c/Users/admin/Documents/Projets/Phylogeny_results]

REPLICON=$1
WORK_DIR=$2
THREADS=5
GENOMES_DIR=$WORK_DIR/genomes
PROKKA_FILES=$WORK_DIR/phylogeny/$REPLICON/sequences/*/prokka/*.gff
PANAROO_DIR=$WORK_DIR/phylogeny/$REPLICON/panaroo
SNIPPY_DIR=$WORK_DIR/phylogeny/$REPLICON/snippy
IQTREE_DIR=$WORK_DIR/phylogeny/$REPLICON/iqtree

# WORKDIR |
#         > genomes
#         > phylogeny |
#                     > paeruginosa |
#                                   > sequences 
#                                   > panaroo
#                                   > iqtree
#                                   > snippy

#Step 0- Pan-genome analysis with Panaroo
mkdir -p $PANAROO_DIR
panaroo -i $PROKKA_FILES -o $PANAROO_DIR --clean-mode strict -a core --core_threshold 1 -t $THREADS &> $PANAROO_DIR/panaroo.log

#Step 1- FASTA of Core genome consensus
# Need python script "core_genome.py" to be in the same folder as this script
# Output : FASTA "core_genome_reference.fa" in panaroo results directory
mkdir -p $PANAROO_DIR # normally this folder already exists
python3 core_genome.py -d $PANAROO_DIR

#Step 2- Snippy input (with Snippy-multi)
# Need python script "snippy_multi_list.py" to be in the same folder as this script
SAMPLES=$(basename -a $GENOMES_DIR/* | paste -d ' ' -s -)
echo $SAMPLES
mkdir -p $SNIPPY_DIR
#Creation of input.tab containing path to all samples R1/R2 Illumina Fastq files
python3 snippy_multi_list.py -d $WORK_DIR -l $SAMPLES -o $SNIPPY_DIR

INPUT_TAB=$SNIPPY_DIR/input.tab
REF_FASTA=$PANAROO_DIR/core_genome_reference.fa
#Creation of snippy commands
snippy-multi $INPUT_TAB --ref $REF_FASTA --cpus $THREADS --force > $SNIPPY_DIR/snippy_commands.sh
#Adding absolute path for output files
sed -i -e "s|snippy-core --ref '|snippy-core --prefix $SNIPPY_DIR/core --ref '|g" $SNIPPY_DIR/snippy_commands.sh

#Step 3- Running Snippy
sh $SNIPPY_DIR/snippy_commands.sh 1> $SNIPPY_DIR/snippy.log 2> $SNIPPY_DIR/snippy.err

#Step4 - IQTREE
CORE_SNPS_ALN=$SNIPPY_DIR/core.full.aln
mkdir -p $IQTREE_DIR
iqtree -s $CORE_SNPS_ALN -m HKY --prefix $IQTREE_DIR/$REPLICON -T $THREADS 1> $IQTREE_DIR/iqtree.log 2> $IQTREE_DIR/iqtree.err

