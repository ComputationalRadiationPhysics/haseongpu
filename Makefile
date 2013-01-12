#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"


# compiler, linker, archiver
NVCC = nvcc

# build variables
SRCS = $(wildcard src/*.cu src/*/*.cu)


all: octrace
	notify-send "COMPILATION" "finished"
	

# build octrace
octrace:  
	$(NVCC) $(SRCS) -o bin/octrace	

