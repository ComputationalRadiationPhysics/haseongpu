#! /usr/bin/env sh
#PBS -q k20
#PBS -l nodes=1:ppn=2
#PBS -l walltime=00:30:00
#PBS -N calcPhiASE_postprocess

cd ~/octrace

rm -f postprocess
touch postprocess
