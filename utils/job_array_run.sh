#! /usr/bin/env sh


## PARAMETER ########################################
MIN_SAMPLE_I=0
MAX_SAMPLE_I=$(($1 - 1))
SUM=0
PIPE_FINISHED="tmp/octrace_job_array_pipe_finished"
PIPE_STARTED="tmp/octrace_job_array_pipe_started"
HOSTNAMES="tmp/hostnames"
STDOUT_PATH="tmp/"
POSTPROCESS="job_array_post.sh"
SUBMIT="utils/submit_k20.sh"
#####################################################

echo " "
cd ~/octrace

echo "Prepare environment..."
echo " "
mkdir -p $STDOUT_PATH
rm -f $PIPE_FINISHED
rm -f $PIPE_STARTED
rm -f $HOSTNAMES
touch $PIPE_STARTED
touch $PIPE_FINISHED
touch $HOSTNAMES
echo 0 >> $PIPE_FINISHED
echo 0 >> $PIPE_STARTED
echo `hostname` >> $HOSTNAMES

echo "Build..."
make
echo " "
echo "minSample: $MIN_SAMPLE_I"
echo "maxSample: $MAX_SAMPLE_I"

echo " "
echo "Submit jobs..."
JOBNAME=`qsub -t $MIN_SAMPLE_I-$MAX_SAMPLE_I $SUBMIT -e $STDOUT_PATH -o $STDOUT_PATH`
#qsub -W depend=afterok:$JOBNAME $POSTPROCESS -e $STDOUT_PATH -o $STDOUT_PATH

echo " "
echo "Wait for jobs..."

MAX_SAMPLE_I=$(($MAX_SAMPLE_I + 1))
while [ $SUM -lt $MAX_SAMPLE_I ] ; do 
    SUM=0;while read l; do SUM=$((SUM+$l));done<$PIPE_FINISHED;
    SUM_START=0;while read l; do SUM_START=$((SUM_START+$l));done<$PIPE_STARTED;
    printf "%d of %d samples finished(%d started)\r" $SUM $MAX_SAMPLE_I $SUM_START
    sleep 1s
done;


rm $PIPE_FINISHED
rm $PIPE_STARTED
echo " "
echo "Done..."