#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"
export CPLUS_INCLUDE_PATH=

# compiler, linker, archiver
NVCC = nvcc

# build flags
LIBS =  -lpthread -lcudart -lm
ARCH = -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_35,code=sm_35
 #NVCC_FLAGS = --use_fast_math -Xptxas="-v"
 #DEBUG_FLAGS = -g -G -lineinfo -D THRUST_DEBUG
NVCC_FLAGS = --use_fast_math "-DBOOST_NOINLINE=__attribute__((noinline))"
GCC_FLAGS = -std=c++0x -J 8 -O2
#DEV_FLAGS = --compiler-options="-Wall -Wextra" 
DEV_FLAGS = --compiler-options="-Wextra -Waddress -Warray-bounds -Wchar-subscripts -Wcomment -Wformat -Wmain -Wmissing-braces -Wparentheses -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-overflow=1 -Wswitch -Wtrigraphs -Wuninitialized -Wunused-function -Wunused-label -Wunused-value -Wunused-variable -Wvolatile-register-var -Werror" 
#FINAL_BUILD = -D NDEBUG


# build variables
SRCS = $(wildcard src/*.cu src/*/*.cu)
OBJS = $(SRCS:src/%.cu=bin/%.o)
INCLUDES = include

all:	
	@mkdir -p bin
	@mkdir -p example/matlab_example/bin
	@mkdir -p example/c_example/bin
	@mkdir -p example/c_example/output
	@cp src/calcPhiASE.m .
	@cp src/calcPhiASE.m example/matlab_example/
	@make calcPhiASE

bin/calc_phi_ase_mpi.o: src/calc_phi_ase_mpi.cc include/calc_phi_ase_mpi.hpp
	@CPLUS_INCLUDE_PATH=/opt/pkg/devel/cuda/5.0/include mpic++ -Wall -Wextra -Werror -lm -c $< -I $(INCLUDES) -o bin/calc_phi_ase_mpi.o

bin/%.o: src/%.cu $(wildcard include/*.hpp)
	@$(NVCC) -dc $< -odir bin --include-path $(INCLUDES) $(ARCH) $(NVCC_FLAGS) $(DEV_FLAGS)  $(DEBUG_FLAGS) $(FINAL_BUILD)

calcPhiASE: $(OBJS) Makefile bin/calc_phi_ase_mpi.o
	@rm -f bin/link.o
	@$(NVCC) $(ARCH) bin/*.o -dlink -o bin/link.o
	@mpic++ bin/*.o -o bin/calcPhiASE $(GCC_FLAGS) $(LIBS)
	@cp bin/calcPhiASE example/matlab_example/bin
	@cp bin/calcPhiASE example/c_example/bin

clean:
	@rm -f bin/*

new: 
	@make clean
	@make

archive:
	rm -f program_sources.zip
	zip -r program_sources.zip src
	zip -r program_sources.zip example
	zip -r program_sources.zip CMakeLists.txt
	zip -r program_sources.zip include
	zip -r program_sources.zip LICENSE.md
	zip -r program_sources.zip README.md
	zip -r program_sources.zip REFERENCE.md
	zip -r program_sources.zip COPYING
	zip -r program_sources.zip utils
