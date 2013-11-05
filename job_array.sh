#! /usr/bin/env sh

MIN_SAMPLE_I=0
MAX_SAMPLE_I=$(($1 - 1))
SUM=0
PIPE_FINISHED="tmp/octrace_job_array_pipe_finished"
PIPE_STARTED="tmp/octrace_job_array_pipe_started"
STDOUT_PATH="tmp/"

echo " "
cd ~/octrace

echo "Prepare environment..."
echo " "
mkdir -p $STDOUT_PATH
rm -f $PIPE_FINISHED
touch $PIPE_FINISHED
rm -f $PIPE_STARTED
touch $PIPE_STARTED
echo 0 >> $PIPE_FINISHED
echo 0 >> $PIPE_STARTED

echo "Build..."
make
echo " "
echo "minSample: $MIN_SAMPLE_I"
echo "maxSample: $MAX_SAMPLE_I"

echo " "
echo "Submit jobs..."
JOBNAME=`qsub -V -t $MIN_SAMPLE_I-$MAX_SAMPLE_I utils/submit_k20.sh -e $STDOUT_PATH -o $STDOUT_PATH`

echo " "
echo "Wait for jobs..."
#qsub analyze.sh -W depend=$JOBNAME:427[]

MAX_SAMPLE_I=$(($MAX_SAMPLE_I + 1))
while [ $SUM -lt $MAX_SAMPLE_I ] ; do 
    SUM=0;while read l; do SUM=$((SUM+$l));done<$PIPE_FINISHED;
    SUM_START=0;while read l; do SUM_START=$((SUM_START+$l));done<$PIPE_STARTED;
    printf "%d of %d samples finished(%d started)\r" $SUM $MAX_SAMPLE_I SUM_START
    sleep 1s
done;


rm $PIPE_FINISHED
rm $PIPE_STARTED
echo " "
echo "Done..."