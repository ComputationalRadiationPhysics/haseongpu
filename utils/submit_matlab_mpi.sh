#!/bin/bash
#PBS -q k20
#PBS -l nodes=10:ppn=8
#PBS -l walltime=96:00:00
#PBS -N calcPhiASE
#PBS -d .
#PBS -o output
#PBS -e output

## ENVIRONMENT ############################
. /opt/modules-3.2.6/Modules/3.2.6/init/bash
export MODULES_NO_OUTPUT=1
module load ~/own.modules
module load analysis/matlab
export -n MODULES_NO_OUTPUT

cd ~/octrace/utils/benchmark_refl

RANKFILE=/tmp/ase_flux_mpirank
mpiexec -npernode 4 --output-filename $RANKFILE echo 1

if [ -e "${RANKFILE}.1.0" ] || [ -e "${RANKFILE}.1.00" ] ; then
  echo found headnode on Host 
  hostname
  rm -f ${RANKFILE}*
  matlab -r "global Runmode; Runmode='mpi'; ASE_calc";
else
  rm -f ${RANKFILE}*
fi
