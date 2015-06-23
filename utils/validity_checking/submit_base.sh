#!/bin/bash
#PBS -q k20
#PBS -l nodes=2:ppn=8
#PBS -l walltime=12:00:00
#PBS -N HASE_integration_validity

## ENVIRONMENT ############################
. /opt/modules-3.2.6/Modules/3.2.6/init/bash
export MODULES_NO_OUTPUT=1
module load ~/own.modules.kepler
export -n MODULES_NO_OUTPUT
uname -a
echo " "
