#!/bin/bash
#PBS -q k20
#PBS -l nodes=1:ppn=8
#PBS -l walltime=20:30:00
#PBS -N calcPhiASE
#PBS -d .
#PBS -o output/
#PBS -e output/

## ENVIRONMENT ############################
. /opt/modules-3.2.6/Modules/3.2.6/init/bash
export MODULES_NO_OUTPUT=1
module load ~/own.modules
module load analysis/matlab
export -n MODULES_NO_OUTPUT

cd ~/octrace/utils/benchmark_refl
matlab -r ASE_calc
