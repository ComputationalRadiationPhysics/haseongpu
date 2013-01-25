#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"


# compiler, linker, archiver
NVCC = nvcc
NVCC_FLAGS = --ptxas-options=-v --use_fast_math --include-path include
ARCH = -arch=sm_30
# --maxrregcount=40

# build variables
SRCS = $(wildcard src/*.cu src/*/*.cu)
INCLUDES = include

all: octrace

# build octrace
octrace:
	$(NVCC) $(SRCS) -dc --include-path $(INCLUDES) $(ARCH) $(NVCC_FLAGS)
	$(NVCC) $(ARCH) *.o -dlink -o link.o
	g++ *.o -o bin/octrace -lcudart
	rm *.o
