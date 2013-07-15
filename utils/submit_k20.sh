#!/bin/bash
#PBS -q k20f
#PBS -l nodes=1:ppn=2
#PBS -l walltime=3:30:00
#PBS -N octrace

. /opt/modules-3.2.6/Modules/3.2.6/init/bash
export MODULES_NO_OUTPUT=1
module load ~/own.modules
export -n MODULES_NO_OUTPUT

uname -a

echo
cd ~/octrace

make

MODE="for_loops"
RAYS=1000000
echo "RAYS: $RAYS"
echo "MODE: $MODE"

FOLDER="$(pwd)"
echo "Executing..."
echo
time ./bin/octrace --mode=$MODE --rays=$RAYS --experiment="$FOLDER/utils/testdata_2" --silent
