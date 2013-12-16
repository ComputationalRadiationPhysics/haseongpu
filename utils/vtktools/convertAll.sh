#!/bin/bash

for i in *.vtk ; do
  convertVTKtoVTU $i
done


