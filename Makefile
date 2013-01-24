#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"


# compiler, linker, archiver
NVCC = nvcc
NVCC_FLAGS = -arch=sm_30 --use_fast_math --include-path include
# --maxrregcount=40

# build variables
SRCS = $(wildcard src/*.c src/*.cu src/*/*.cu)
INCLUDES = include

all: octrace

# build octrace
octrace:  
	$(NVCC) $(SRCS) -dc --include-path $(INCLUDES) $(NVCC_FLAGS)
	$(NVCC) -arch=sm_30 *.o -dlink -o link.o
	g++ *.o -o bin/octrace -lcudart
	rm *.o
