# Launch parallel tasks on a list of samples
# Give task name then sample list (in current directory)
# eg, launch.sh flyepilon sample.txt

while read SAMPLE || [ -n "$SAMPLE" ]; do
  /bin/bash "$1".sh $SAMPLE &
done < ./$2

# Snippet for batch renaming if needed
# for i in ./*/amrfinder/amrfinder.err; do mv $i $(sed 's/err/log/g' <<<$i); done
# remark sed <<<$i idiom to avoid reading the file pointed by $i

# SOLVED 2022-10-24
# Last sample in text file not executed, possibly due to missing line end ?
# See https://stackoverflow.com/questions/12916352/shell-script-read-missing-last-line
# Solution is the || [ -n "$SAMPLE" ] block to prevent read to throw error
# at EOF and stop the loop.