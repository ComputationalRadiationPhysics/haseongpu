#! /usr/bin/env sh

if [ "$1" = "" ] || [ "$2" = "" ]
    then
    echo "Usage: $0 NUM_SAMPLES NUM_NODES"
    exit 1
fi


## PARAMETER ########################################
MIN_SAMPLE_I=0
MAX_SAMPLE_I=$(($1 - 1))
NUM_SAMPLES=$1
NUM_NODES=$2
SUM=0
TMP_PATH="tmp"
PIPE_FINISHED="$TMP_PATH/octrace_job_array_pipe_finished"
PIPE_STARTED="$TMP_PATH/octrace_job_array_pipe_started"
HOSTNAMES="$TMP_PATH/hostnames"
RESULTS="output/results/"
POSTPROCESS="job_array_post.sh"
SUBMIT="utils/submit_k20.sh"
#####################################################

echo " "
cd ~/octrace

echo "Prepare environment..."
echo " "
rm -f $TMP_PATH/*
mkdir -p $RESULTS
mkdir -p $TMP_PATH
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
JOBNAME=`qsub -t 0-$(($NUM_NODES-1)) $SUBMIT -e $TMP_PATH/ -o $TMP_PATH/ -v NUM_NODES=$NUM_NODES,NUM_SAMPLES=$NUM_SAMPLES`
#qsub -W depend=afterok:$JOBNAME $POSTPROCESS -e $TMP_PATH/ -o $TMP_PATH/

echo " "
echo "Wait for jobs..."

MAX_SAMPLE_I=$(($MAX_SAMPLE_I + 1))
while [ $SUM -lt $NUM_NODES ] ; do 
    SUM=0;while read l; do SUM=$((SUM+$l));done<$PIPE_FINISHED;
    SUM_START=0;while read l; do SUM_START=$((SUM_START+$l));done<$PIPE_STARTED;
    printf "%d of %d nodes finished(%d started)\r" $SUM $NUM_NODES $SUM_START
    sleep 1s
done;

## CAT RESULTS ################################
rm -f $TMP_PATH/phi_ASE.txt
for i in $(ls $RESULTS)
do
    cat $RESULTS/$i >> output/phi_ASE.txt
done


rm $PIPE_FINISHED
rm $PIPE_STARTED
echo " "
echo "Done..."