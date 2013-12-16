#!/bin/bash

if [ "$1" = "" ] ; then
  echo "give a pattern for the VTK-files! (like dndt_ASE, if the files are called dndt_ASE_XX.vtk)"
  exit 1
fi

PATTERN="$1"
FILE="animation_$1.pvd"

echo '<?xml version"1.0"?>' > $FILE
echo '<VTKFile type="Collection" version="0.1">' >> $FILE
echo '  <Collection>' >> $FILE

j=1;


for i in $(ls -v $PATTERN*.vtk)
do
  echo converting $i
  ./convertVTKtoVTU $i
  OUTPUT="$(basename -s .vtk $i)"
  echo '    <DataSet timestep="'$j'" file="'$OUTPUT'.vtu"/>' >> $FILE
  j=$(($j+1))
done

echo '  </Collection>' >> $FILE
echo '</VTKFile>' >> $FILE
