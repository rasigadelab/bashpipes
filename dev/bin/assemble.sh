# Launch parallel tasks for assembling

while read SAMPLE; do
  flyepilon.sh $SAMPLE &
done < ./$1