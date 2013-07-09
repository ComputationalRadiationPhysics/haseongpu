#!/bin/bash

MODE="2"
#MODE="ray_propagation_gpu"
EXPERIMENT="testdata_2"
SILENT="--silent"
RAYSPERSAMPLE="10000"
COPY_TO_REMOTE=0
REMOTE_USER="s0500037"
REMOTE_HOST="ganymed.inf.tu-dresden.de"
REMOTE_FOLDER="/home/$REMOTE_USER/."

MODE="for_loops"
#MODE="ray_propagation_gpu"

cd ~/octrace
FOLDER="$(pwd)"
make && ./bin/octrace --experiment="$FOLDER/utils/$EXPERIMENT" --mode=$MODE $SILENT --rays=$RAYSPERSAMPLE
#./bin/octrace --experiment=~/octrace/utils/testdata_2/ --silent --rays=10000 --mode=2

if [ $COPY_TO_REMOTE -eq 1 ] ; then
	echo "BASH: Copy octrace.vtk to remote host"
	scp octrace.vtk $REMOTE_USER@$REMOTE_HOST:$REMOTE_FOLDER
fi
