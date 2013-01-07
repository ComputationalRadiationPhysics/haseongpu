#
#	Makefile
#
DATE="`date +%y%m%d%H%M%S`"

# commandline
CMD = -O0 -g -J2
ARGS = 
SPACE = " "

# compiler, linker, archiver
CPP = g++
#CPP = colorgcc
#CPP = clang
DOXYGEN = doxygen

# compiler flags
LIBS		= 
CPPINCLUDES 	= -I./include 
COMMON_CPPFLAGS = $(CPPINCLUDES)
CPPFLAGS 	= $(COMMON_CPPFLAGS) -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -g3  -fno-strict-aliasing -g 
LDFLAGS 	= -L. 


# build variables
SRCS = $(wildcard src/*.cc src/*/*.cc)
OBJS = $(SRCS:.cc=.o)
DEPS = $(SRCS:.cc=.d)

# for building the component libraries
DEPS += $(LIBSRCS:.cc=.d)

all: octrace
	notify-send "COMPILATION" "finished"	

# build octrace
octrace: $(OBJS)
	$(CPP) -o $@ $(OBJS) $(LIBS) $(LDFLAGS) $(CMD) $(CPPFLAGS) $(ARGS)

# build object file and dependencies files
.cc.o:
	$(CPP) $(CMD) -MD $(CPPFLAGS) $(ARGS) -c -o $@ $<

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