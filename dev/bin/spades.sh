# Illumina de novo assembly

# SPADES ASSEMBLER

SAMPLE=$1
SAMPLE_R1="$SAMPLE""_R1.fastq.gz"
SAMPLE_R2="$SAMPLE""_R2.fastq.gz"
SAMPLE_POLISHED="$SAMPLE""_polished"
SAMPLE_POLISHED_REALIGNED="$SAMPLE""_polished_realigned"

cd $SAMPLE
mkdir spades
mkdir quast
mkdir circlator

spades.py -1 $SAMPLE_R1 -2 $SAMPLE_R2 --isolate -o spades -t 48 \
    1> spades/spades_1.log 2> spades/spades.err && \
cp spades/contigs.fasta ./"$SAMPLE_POLISHED".fasta && \
quast.py -o quast ./"$SAMPLE_POLISHED".fasta \
    1> quast/quast.log \
    2> quast/quast.err && \
circlator fixstart "$SAMPLE_POLISHED".fasta circlator/"$SAMPLE_POLISHED_REALIGNED" && \
cp circlator/$SAMPLE_POLISHED_REALIGNED.fasta ./"$SAMPLE_POLISHED_REALIGNED".fasta