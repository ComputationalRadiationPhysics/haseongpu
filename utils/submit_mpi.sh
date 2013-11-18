#PBS -q k20
#PBS -l nodes=1:ppn=8
#PBS -N ase_flux_MPI
#PBS -l walltime=01:00:00
#PBS -d .

# SETTINGS #################
MAXGPUS="1"
RAYSPERSAMPLE="10000"
MAXRAYS="10000000"
MIN_SAMPLE_I=0
MAX_SAMPLE_I=3209
MODE="mpi"
GPU_PER_NODE="2"
############################

cd ~/octrace
FOLDER="$(pwd)"

mpiexec -npernode $GPU_PER_NODE ./bin/calcPhiASE --experiment="$FOLDER/input" --mode=$MODE $SILENT --rays=$RAYSPERSAMPLE --compare="$FOLDER/$COMPARE" $WRITE_VTK $USE_REFLECTION --maxrays=$MAXRAYS --maxgpus=$MAXGPUS --min_sample_i=$MIN_SAMPLE_I --max_sample_i=$MAX_SAMPLE_I 
