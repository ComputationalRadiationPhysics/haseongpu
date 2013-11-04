#! /usr/bin/env sh

MIN_SAMPLE_I=0
MAX_SAMPLE_I=$(($1 - 1))
SUM=0
PIPE="job_array_pipe"
STDOUT_PATH="tmp/"

echo " "
cd ~/octrace

echo "Prepare environment..."
echo " "
mkdir -p $STDOUT_PATH
rm -f $PIPE
touch $PIPE
echo 0 >> $PIPE

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
    SUM=0;while read l; do SUM=$((SUM+$l));done<$PIPE;
    printf "%d of %d samples finished\r" $SUM $MAX_SAMPLE_I
    sleep 1s
done;


rm $PIPE
echo " "
echo "Done..."