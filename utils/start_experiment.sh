#!/bin/bash

MODE="2"
#MODE="ray_propagation_gpu"
EXPERIMENT="testdata_2"
SILENT="--silent"
RAYSPERSAMPLE="10000"
COPY_TO_REMOTE=1
REMOTE_USER="s0500037"
REMOTE_HOST="ganymed.inf.tu-dresden.de"
REMOTE_FOLDER="/home/$REMOTE_USER/."

MODE="ray_propagation_gpu"
REFLECTION="--reflection"
MAXRAYSPERSAMPLE="10000000"
EXPECTATIONTHRESHOLD="0.005"

#COMPARE="octrace_backtotheroots_100M.vtk"

cd ~/octrace
FOLDER="$(pwd)"
make && ./bin/octrace --experiment="$FOLDER/utils/$EXPERIMENT" --mode=$MODE $SILENT --rays=$RAYSPERSAMPLE --compare=$COMPARE $REFLECTION --expectation=$EXPECTATIONTHRESHOLD --maxrays=$MAXRAYSPERSAMPLE --write-vtk

if [ $COPY_TO_REMOTE -eq 1 ] ; then
	echo "BASH: Copy octrace.vtk to remote host"
	scp octrace.vtk $REMOTE_USER@$REMOTE_HOST:$REMOTE_FOLDER
fi
