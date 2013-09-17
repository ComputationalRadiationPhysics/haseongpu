#!/bin/bash
#PBS -q k20
#PBS -l nodes=1:ppn=2
#PBS -l walltime=40:30:00
#PBS -N octrace_baseline_upto_100M_rays_per_sample_14r

. /opt/modules-3.2.6/Modules/3.2.6/init/bash
export MODULES_NO_OUTPUT=1
module load ~/own.modules
export -n MODULES_NO_OUTPUT

uname -a

echo
cd ~/octrace

make

EXPERIMENT="testdata_2"
SILENT="--silent"
DEVICE=0
#USE_REFLECTION="--reflection"
EXPECTATION="0.001"
MAXRAYS="100000000"
WRITE_VTK="--write-vtk"
RAYSPERSAMPLE="1000000"
MODE="ray_propagation_gpu"
echo "RAYS: $RAYSPERSAMPLE"
echo "MODE: $MODE"


FOLDER="$(pwd)"
echo "Executing..."
echo
time ./bin/octrace --experiment="$FOLDER/utils/$EXPERIMENT" --mode=$MODE $SILENT --rays=$RAYSPERSAMPLE --compare="$FOLDER/$COMPARE" --device="$DEVICE" $WRITE_VTK $USE_REFLECTION --expectation=$EXPECTATION --maxrays=$MAXRAYS
