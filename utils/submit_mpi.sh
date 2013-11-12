#PBS -q k20
#PBS -l nodes=4:ppn=8
#PBS -N ase_flux_MPI
#PBS -l walltime=00:10:00
#PBS -d .

cd ~/octrace

mpiexec -npernode 4 ~/octrace/bin/mpi_ase
