#!/bin/bash

EXPERIMENT="testdata_2"
SILENT="--silent"
RAYSPERSAMPLE="1000000"
COPY_TO_REMOTE=1
REMOTE_USER="s0500037"
REMOTE_HOST="ganymed.inf.tu-dresden.de"
REMOTE_FOLDER="/home/$REMOTE_USER/."

cd ~/octrace
FOLDER="$(pwd)"
make && ./bin/octrace --experiment="$FOLDER/utils/$EXPERIMENT" --mode=ray_propagation_gpu $SILENT --rays=$RAYSPERSAMPLE

if [ $COPY_TO_REMOTE -eq 1 ] ; then
	echo "BASH: Copy octrace.vtk to remote host"
	scp octrace.vtk $REMOTE_USER@$REMOTE_HOST:$REMOTE_FOLDER
fi
