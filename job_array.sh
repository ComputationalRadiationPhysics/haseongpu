#! /usr/bin/env sh

MIN_SAMPLE_I=0
MAX_SAMPLE_I=$(($1 - 1))
FIFO="job_array_fifo"
SUM=0

echo " "
cd ~/octrace

echo "Build..."
make
echo " "
echo "minSample: $MIN_SAMPLE_I"
echo "maxSample: $MAX_SAMPLE_I"

rm -f $FIFO
mkfifo $FIFO

qsub -V -t $MIN_SAMPLE_I-$MAX_SAMPLE_I utils/submit_k20.sh 

#while [ $SUM -lt $MAX_SAMPLE_I ] ; do 
#    TMP=`cat $FIFO`
#    SUM=$(($SUM + $TMP))
#    echo $SUM
#done;

echo " "
echo "Done..."