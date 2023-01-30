# Illumina de novo assembly

# SPADES ASSEMBLER

SAMPLE=$1
SAMPLE_R1="$SAMPLE""_R1.fastq.gz"
SAMPLE_R2="$SAMPLE""_R2.fastq.gz"
SAMPLE_POLISHED="$SAMPLE""_polished"

cd $SAMPLE
mkdir spades
mkdir quast

spades.py -1 $SAMPLE_R1 -2 $SAMPLE_R2 --isolate -o spades \
    1> spades/spades_1.log 2> spades/spades.err && \
cp spades/contigs.fasta $SAMPLE_POLISHED && \
quast.py -o quast $SAMPLE_POLISHED \
    1> quast/quast.log \
    2> quast/quast.err 