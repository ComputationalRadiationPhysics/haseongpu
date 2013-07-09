#!/bin/bash

cd ~/octrace/bin

for i in {4,5,6,7} 
do
  R=$(( 10 ** $i ))
  qsub -v RAYS="$R",MODE="ray_propagation_gpu" ~/octrace/utils/submit_k20.sh
done


for i in {4,5} 
do
  R=$(( 10 ** $i ))
  qsub -v RAYS="$R",MODE="for_loops" ~/octrace/utils/submit_k20.sh
done
