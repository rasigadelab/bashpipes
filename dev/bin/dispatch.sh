while read SAMPLE; do
    mkdir $SAMPLE &&
    mv "$SAMPLE"*.fastq.gz ./"$SAMPLE"
done < ./samples.txt