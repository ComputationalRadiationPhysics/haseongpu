#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"


# compiler, linker, archiver
NVCC = nvcc
NVCC_FLAGS = --use_fast_math
DEV_FLAGS = --compiler-options="-Wall -Wextra"
#NVCC_FLAGS = --use_fast_math --include-path include
ARCH = -arch=sm_20
#ARCH = -arch=sm_35

# --maxrregcount=40

# build variables
SRCS = $(wildcard src/*.cu src/*/*.cu)
OBJS = $(SRCS:src/%.cu=bin/%.o)
TESTSRCS = $(wildcard tests/*.cu)
TEST_FLAGS = -g -G
INCLUDES = include

all: octrace

bin/%.o: src/%.cu $(wildcard include/*.h)
	$(NVCC) -dc $< -odir bin --include-path $(INCLUDES) $(ARCH) $(NVCC_FLAGS) $(DEV_FLAGS)

octrace: $(OBJS) Makefile
	rm -f bin/link.o
	$(NVCC) $(ARCH) bin/*.o -dlink -o bin/link.o
	g++ bin/*.o -o bin/octrace -lcudart

clean:
	rm -f bin/*

new: 
	make clean
	make

final_build:
	rm -f bin/link.o
	$(NVCC) $(SRCS) -dc -odir bin --include-path $(INCLUDES) $(ARCH) $(NVCC_FLAGS)
	$(NVCC) $(ARCH) bin/*.o -dlink -o bin/link.o
