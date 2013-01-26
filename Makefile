#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"


# compiler, linker, archiver
NVCC = nvcc
#NVCC_FLAGS = --ptxas-options=-v --use_fast_math --include-path include
NVCC_FLAGS = --use_fast_math --include-path include
ARCH = -arch=sm_30
# --maxrregcount=40

# build variables
SRCS = $(wildcard src/*.cu src/*/*.cu)
SRCS = $(wildcard tests/*.cu)
INCLUDES = include

all: octrace

# build octrace
octrace:
	$(NVCC) $(SRCS) -dc -odir bin --include-path $(INCLUDES) $(ARCH) $(NVCC_FLAGS)
	$(NVCC) $(ARCH) bin/*.o -dlink -o bin/link.o
	g++ bin/*.o -o bin/octrace -lcudart

test: 
	$(NVCC) $(SRCS) $(TESTSRCS) -dc -odir bin --include-path $(INCLUDES) $(ARCH) $(NVCC_FLAGS)
	rm bin/main.o
	$(NVCC) $(ARCH) bin/*.o -dlink -o bin/link.o
	g++ bin/*.o -o bin/octrace -lcudart
