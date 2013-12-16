#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"
export CPLUS_INCLUDE_PATH=

# compiler, linker, archiver
NVCC = nvcc
NVCC_FLAGS = --use_fast_math -Xptxas="-v"
NVCC_FLAGS = --use_fast_math
GCC_FLAGS = -std=c++0x -J 8 -O2
LIBS =  -lpthread -lcudart -lm

DEV_FLAGS = --compiler-options="-Wall -Wextra"
THRUST_FLAGS = -D THRUST_DEBUG

ARCH = -arch=sm_20
ARCH = -arch=sm_35
ARCH = -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_35,code=sm_35


# --maxrregcount=40

# build variables
SRCS = $(wildcard src/*.cu src/*/*.cu)
OBJS = $(SRCS:src/%.cu=bin/%.o)
TESTSRCS = $(wildcard tests/*.cu)
DEBUG_FLAGS = -g -G -lineinfo
INCLUDES = include

all: calcPhiASE

bin/calc_phi_ase_mpi.o: src/calc_phi_ase_mpi.cc include/calc_phi_ase_mpi.h
	CPLUS_INCLUDE_PATH=/opt/pkg/devel/cuda/5.0/include mpic++ -Wall -Wextra -lm -c $< -I include -o bin/calc_phi_ase_mpi.o

bin/%.o: src/%.cu $(wildcard include/*.h)
	$(NVCC) -dc $< -odir bin --include-path $(INCLUDES) $(ARCH) $(NVCC_FLAGS) $(DEV_FLAGS) 
#$(DEBUG_FLAGS) $(THRUST_FLAGS)

calcPhiASE: $(OBJS) Makefile bin/calc_phi_ase_mpi.o
	rm -f bin/link.o
	mkdir -p bin
	mkdir -p output
	mkdir -p output/calcPhiAseTmp
	mkdir -p output/vtk
	mkdir -p input
	$(NVCC) $(ARCH) bin/*.o -dlink -o bin/link.o
	mpic++ bin/*.o -o bin/calcPhiASE $(GCC_FLAGS) $(LIBS)
	cp src/calcPhiASE.m .

clean:
	rm -f bin/*

new: 
	make clean
	make

final_build:
	rm -f bin/link.o
	mkdir -p bin
	$(NVCC) $(SRCS) -dc -odir bin --include-path $(INCLUDES) $(ARCH) $(NVCC_FLAGS)
	$(NVCC) $(ARCH) bin/*.o -dlink -o bin/link.o
	cp src/calcPhiASE.m .

#mpi: src/calc_phi_ase_mpi.cc include/calc_phi_ase_mpi.h
#	CPLUS_INCLUDE_PATH=/opt/pkg/devel/cuda/5.0/include mpic++ -Wall -Wextra -lm -c $< -I include -o bin/calc_phi_ase_mpi.o

