#! /usr/bin/env sh
#PBS -q k20
#PBS -l nodes=1:ppn=2
#PBS -l walltime=00:30:00
#PBS -N calcPhiASE

# INIT ENVIRONMENT ###############################
. /opt/modules-3.2.6/Modules/3.2.6/init/bash
export MODULES_NO_OUTPUT=1
#module load ~/own.modules
export -n MODULES_NO_OUTPUT
uname -a

echo
cd ~/octrace

## OCTRACE PARAMETER #############################
MAXGPUS="1"
#USE_REFLECTION="--reflection"
EXPECTATION="0.005"
RAYSPERSAMPLE="10000"
MAXRAYS="10000"
EXPERIMENT="testdata_2"
SILENT="--silent"
SAMPLE=$PBS_ARRAYID
MODE="ray_propagation_gpu"

## CLUSTER PARAMETER #############################
PIPE_FINISHED="tmp/octrace_job_array_pipe_finished"
PIPE_STARTED="tmp/octrace_job_array_pipe_started"
HOSTNAMES="tmp/hostnames"

FOLDER="$(pwd)"
echo "Executing..."
echo

echo `hostname` >> $HOSTNAMES
echo 1 >> $PIPE_STARTED

time ./bin/calcPhiASE --experiment="$FOLDER/utils/$EXPERIMENT" --mode=$MODE $SILENT --rays=$RAYSPERSAMPLE $USE_REFLECTION --maxrays=$MAXRAYS --maxgpus=$MAXGPUS --sample_i=$SAMPLE

echo 1 >> $PIPE_FINISHED
