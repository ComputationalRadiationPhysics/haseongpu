#!/bin/bash
#PBS -q k20
#PBS -l nodes=1:ppn=8
#PBS -l walltime=00:30:00
#PBS -N calcPhiASE

## ENVIRONMENT ############################
. /opt/modules-3.2.6/Modules/3.2.6/init/bash
export MODULES_NO_OUTPUT=1
module load ~/own.modules
export -n MODULES_NO_OUTPUT
uname -a
echo " "
cd ~/octrace

## FS PARAMETER ###########################
PIPE_STARTED="tmp/octrace_job_array_pipe_started"
PIPE_FINISHED="tmp/octrace_job_array_pipe_finished"
HOSTNAMES="tmp/hostnames"
FOLDER="$(pwd)"
NODE_ID=$PBS_ARRAYID
SAMPLE_PER_NODE="$(echo "$NUM_SAMPLES / $NUM_NODES" | bc -l)"

## OCTRACE PARAMETER ######################
MAXGPUS="4"
#USE_REFLECTION="--reflection"
RAYSPERSAMPLE="100000000"
MAXRAYS="100000000"
SILENT="--silent"

MODE="ray_propagation_gpu"
MIN_SAMPLE_I=$(echo "$NODE_ID * $SAMPLE_PER_NODE" | bc)
MAX_SAMPLE_I=$(echo "(($NODE_ID+1) * $SAMPLE_PER_NODE)/1" | bc)
if [ $MAX_SAMPLE_I -gt $NUM_SAMPLES ]; then MAX_SAMPLE_I=$NUM_SAMPLES; fi
###########################################

echo "Executing..."
echo

echo `hostname` >> $HOSTNAMES
echo 1 >> $PIPE_STARTED

time ./bin/calcPhiASE --experiment="$FOLDER/input" --mode=$MODE $SILENT --rays=$RAYSPERSAMPLE $USE_REFLECTION --maxrays=$MAXRAYS --maxgpus=$MAXGPUS --min_sample_i=$MIN_SAMPLE_I --max_sample_i=$(($MAX_SAMPLE_I-1))

echo 1 >> $PIPE_FINISHED
