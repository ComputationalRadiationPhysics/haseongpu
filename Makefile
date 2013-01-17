#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"


# compiler, linker, archiver
NVCC = nvcc
NVCC_FLAGS = --ptxas-options=-v -arch=sm_30 --use_fast_math --include-path include
# --maxrregcount=40

# build variables
SRCS = $(wildcard src/*.c src/*.cu src/*/*.cu)

all: octrace

# build octrace
octrace:  
	$(NVCC) $(SRCS) -o bin/octrace $(NVCC_FLAGS)

