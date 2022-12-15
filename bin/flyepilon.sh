# Pilon-polished Nanopore de novo assembly

# FLYE ASSEMBLER

SAMPLE=$1
SAMPLE_R1="$SAMPLE""_R1.fastq.gz"
SAMPLE_R2="$SAMPLE""_R2.fastq.gz"
SAMPLE_ONT="$SAMPLE""_ONT.fastq.gz"
SAMPLE_ASM_RAW="$SAMPLE""_assembly_raw"
SAMPLE_BAM="$SAMPLE"".bam"
SAMPLE_BAM_SORTED="$SAMPLE"".sorted.bam"
SAMPLE_POLISHED="$SAMPLE""_polished"

cd $SAMPLE
mkdir flye
mkdir polish

cd flye && \
flye --nano-raw ../"$SAMPLE_ONT" -o . --threads 48 \
    1> flye.log 2> flye.err && \
cd .. && \
cp flye/assembly.fasta ./"$SAMPLE_ASM_RAW".fasta && \
cd polish && \
bowtie2-build ../"$SAMPLE_ASM_RAW".fasta index \
    &> bowtie2.index.log && \
bowtie2 -x index -1 ../"$SAMPLE_R1" -2 ../"$SAMPLE_R2" -p 48 2>> bowtie2.map.log |samtools view -bS - > $SAMPLE_BAM && \
samtools sort $SAMPLE_BAM > $SAMPLE_BAM_SORTED && \
samtools index $SAMPLE_BAM_SORTED &> samtools.index.log && \
pilon -Xmx256G --genome ../"$SAMPLE_ASM_RAW".fasta --bam $SAMPLE_BAM_SORTED --threads 48 --changes --output $SAMPLE_POLISHED \
    &> pilon.log && \
cp $SAMPLE_POLISHED.fasta ../"$SAMPLE_POLISHED".fasta

# https://gist.github.com/stevekm/e054544fe3849bd7173d4c9124577115
# Problem found with samtools sort redirection, steals output from samtools sort => don't log