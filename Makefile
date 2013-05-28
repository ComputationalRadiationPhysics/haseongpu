#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"


# compiler, linker, archiver
NVCC = nvcc
#NVCC_FLAGS = --ptxas-options=-v --use_fast_math --include-path include
NVCC_FLAGS = --use_fast_math --include-path include
ARCH = -arch=sm_35

# --maxrregcount=40

# build variables
SRCS = $(wildcard src/*.cu src/*/*.cu)
OBJS = $(SRCS:src/%.cu=bin/%.o)
TESTSRCS = $(wildcard tests/*.cu)
TEST_FLAGS = -g -G
INCLUDES = include

all: octrace

bin/%.o: src/%.cu $(wildcard include/*.h)
	$(NVCC) -dc $< -odir bin --include-path $(INCLUDES) $(ARCH) $(NVCC_FLAGS)

octrace: $(OBJS) Makefile
	if [ -f bin/link.o ]; then rm bin/link.o; fi;
	$(NVCC) $(ARCH) bin/*.o -dlink -o bin/link.o
	g++ bin/*.o -o bin/octrace -lcudart

#test: 
#	$(NVCC) $(SRCS) $(TESTSRCS) $(TEST_FLAGS) -dc -odir bin --include-path $(INCLUDES) $(ARCH) $(NVCC_FLAGS)
#	rm bin/main.o
#	$(NVCC) $(ARCH) bin/*.o -dlink -o bin/link.o
#	g++ bin/*.o -o bin/test -lcudart
#	rm bin/test*.o
