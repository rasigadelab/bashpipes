#!/bin/bash
# https://stackoverflow.com/questions/356100/how-to-wait-in-bash-for-several-subprocesses-to-finish-and-return-exit-code-0

# Launch a certain number of concurrent tasks, monitor the number of running tasks and launch the next one when done

# TODO next, give explicit names to processes

# TYPICAL EXAMPlE
# & launch.sh 10 dosomestuff samples.txt
# reads arguments in samples.txt, launches script dosomestuff.sh in the background for each argument,
# allowing at most 10 jobs at a time. If the 10-job limit is reached, waits for the first currently running job to finish
# before launching a new one.

# Input arguments
MAX_JOBS=$1
JOB_NAME=$2
QUEUE_FILE=$3

# Log starting information
echo "[$(date +'%0Y-%0m-%0d %0H:%0M:%0S')]" "Starting launcher command:" "$0 $1 $2 $3" >> launch2.log

# SOLVED 2022-10-24
# Last sample in text file not executed, possibly due to missing line end ?
# See https://stackoverflow.com/questions/12916352/shell-script-read-missing-last-line
# Solution is the || [ -n "$SAMPLE" ] block to prevent read to throw error
# at EOF and stop the loop.
while read SAMPLE || [ -n "$SAMPLE" ]; do
    # Launch process and store 
    /bin/bash "$JOB_NAME".sh $SAMPLE &

    # Used RAM in GB
    USED_RAM=$(awk '/^Mem/ {print $3}' <(free -g))
    FREE_RAM=$(awk '/^Mem/ {print $4}' <(free -g))

    # Log command and PID
    echo "[$(date +'%0Y-%0m-%0d %0H:%0M:%0S')]" "Launched PID ""$!"":" "$JOB_NAME".sh $SAMPLE"; used RAM $USED_RAM GB; free RAM $FREE_RAM GB">> launch2.log

    # List all current running processes and make them an array (surrounding parenthesis)
    PIDS=($(jobs -p))

    # Number of currently running jobs
    len=${#PIDS[@]}
    # echo $len

    # Wait for first running job to finish before continuing the loop
    if (( $len >= $MAX_JOBS ))
    then
    # echo "Waiting... for PID "${PIDS[0]}
    echo "[$(date +'%0Y-%0m-%0d %0H:%0M:%0S')]" "Reached job limit, waiting for PID "${PIDS[0]}" to finish." >> launch2.log
    wait ${PIDS[0]}
    fi

done < ./$QUEUE_FILE

echo "[$(date +'%0Y-%0m-%0d %0H:%0M:%0S')]" "All jobs launched, waiting for remaining jobs to finish." >> launch2.log

# Wait for jobs to finish
PIDS=($(jobs -p))
len=${#PIDS[@]}
while (( $len > 0 )); do
    echo "[$(date +'%0Y-%0m-%0d %0H:%0M:%0S')]" "Waiting for PID "${PIDS[0]}" to finish." >> launch2.log
    wait ${PIDS[0]}
    PIDS=($(jobs -p))
    len=${#PIDS[@]}
done


echo "[$(date +'%0Y-%0m-%0d %0H:%0M:%0S')]" "Exiting." >> launch2.log