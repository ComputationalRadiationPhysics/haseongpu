#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"

# commandline
#CMD = -O0 -g -J2 -MD
CMD =
ARGS = 
SPACE = " "

# compiler, linker, archiver
#CPP = g++
#CPP = colorgcc
#CPP = clang
CPP = nvcc
DOXYGEN = doxygen

# compiler flags
LIBS		= 
CPPINCLUDES 	= -I./include 
COMMON_CPPFLAGS = $(CPPINCLUDES)
#CPPFLAGS 	= $(COMMON_CPPFLAGS) -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -g3  -fno-strict-aliasing -g 
CPPFLAGS        =
#LDFLAGS 	= -L. 
LDFLAGS         =


# build variables
SRCS = $(wildcard src/*.cu src/*/*.cu)
OBJS = $(SRCS:.cu=.o)
DEPS = $(SRCS:.cu=.d)

# for building the component libraries
DEPS += $(LIBSRCS:.cu=.d)

all: octrace
	notify-send "COMPILATION" "finished"	

# build octrace
octrace: $(OBJS)
	$(CPP) -o bin/$@ $(OBJS) $(LIBS) $(LDFLAGS) $(CMD) $(CPPFLAGS) $(ARGS)

# build object file and dependencies files
.cu.o:
	$(CPP) $(CMD) $(CPPFLAGS) $(ARGS) -c -o $@ $<

# clean up backups and old files
clean:
	rm -f *~ */*~ */*/*~
	rm -f $(OBJS) $(DEPS) $(LIBOBJS) $(LIBFILES) $(LIBRISSOBJ)
	rm -f log.txt
	rm -f octrace	

# generate documentation
doc:
	$(DOXYGEN) Doxygen.conf

# include headerfile dependencies for sources
-include $(DEPS)
