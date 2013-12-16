#PBS -q k20
#PBS -N ase_flux_MPI_100M
#PBS -l walltime=01:00:00
#PBS -o output/methods_compare
#PBS -e output/methods_compare
#PBS -d .
#PBS -l nodes=14:ppn=8

_10K="10000"
_100K="100000"
_1M="1000000"
_10M="10000000"
_100M="100000000"
_300M="300000000"

cd ~/octrace
FOLDER="$(pwd)"

# SETTINGS #################
MAXGPUS="1"
RAYSPERSAMPLE=${_300M}
MAXRAYS=${_300M}
MIN_SAMPLE_I=0
MAX_SAMPLE_I=3209
MODE="mpi"
GPU_PER_NODE="4"
INPUT="$FOLDER/input"
OUTPUT="$FOLDER/output/"
WRITE_VTK="--write-vtk"
REPETITIONS="1"
############################

mpiexec -npernode $GPU_PER_NODE ./bin/calcPhiASE --input=$INPUT --output=$OUTPUT --mode=$MODE $SILENT --rays=$RAYSPERSAMPLE --compare="$FOLDER/$COMPARE" $WRITE_VTK $USE_REFLECTION --maxrays=$MAXRAYS --maxgpus=$MAXGPUS --min_sample_i=$MIN_SAMPLE_I --max_sample_i=$MAX_SAMPLE_I --repetitions=$REPETITIONS

