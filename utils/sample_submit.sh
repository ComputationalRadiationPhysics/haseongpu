#!/bin/bash
#PBS -q laser
#PBS -l nodes=1:ppn=64
#PBS -l walltime=4:00:00
#PBS -N test

. /opt/modules-3.2.6/Modules/3.2.6/init/bash
export MODULES_NO_OUTPUT=1
module load /home/schulzh/own.modules
export -n MODULES_NO_OUTPUT

uname -a
echo
cd /home/schnei49

echo "Compiling..."
echo
#mpicc -o lasttest rechn.c

echo "Executing..."
echo
echo "Hallo Welt!"
