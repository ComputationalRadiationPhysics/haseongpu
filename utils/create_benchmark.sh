#!/bin/bash

# terminate the script on error
set -e
if [ "$1" = "" ] ; then
  echo usage: 
  echo "$0 <Directory for new benchmark>"
  exit 1
fi

DEST="$1"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


cd $DIR
mkdir $DEST

# create the template structure
cp -r benchmark_template/* $DEST/.

# copy the newest routines to the folder
cp ../src/calcPhiASE.m $DEST/.
cp ../bin/calcPhiASE $DEST/bin/.

# prepare the VTK-Tools
cp vtktools/* $DEST/.
cd $DEST

if [ "$(hostname)" = "hypnos" ] ; then
  make 1>/dev/null
else
  echo "You are not on hypnos, therefore the VTK-tools can not be automatically compiled (that's sad but ok)"
fi
