#!/bin/bash
#PBS -q k20f
#PBS -l nodes=1:ppn=2
#PBS -l walltime=0:10:00
#PBS -N octrace

. /opt/modules-3.2.6/Modules/3.2.6/init/bash
export MODULES_NO_OUTPUT=1
module load ~/own.modules
export -n MODULES_NO_OUTPUT

uname -a

echo
cd ~/octrace

make

echo "Executing..."
echo
./bin/octrace --mode=bruteforce_gpu --rays=256 --experiment=/home/zenker64/octrace/utils/testdata_2
#./bin/octrace --mode=ray_propagation_gpu --rays=256
