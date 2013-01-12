#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"


# compiler, linker, archiver
NVCC = nvcc

# build variables
SRCS = $(wildcard src/*.cu src/*/*.cu)

all: octrace

# build octrace
octrace:  
	$(NVCC) $(SRCS) -o bin/octrace --ptxas-options=-v -arch=sm_30 --use_fast_math

