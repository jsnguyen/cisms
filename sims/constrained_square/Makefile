CC=gcc

_INCDIRS=../../include
INCDIRS=$(addprefix -I,$(_INCDIRS))

_LIBDIRS=../../lib
LIBDIRS=$(addprefix -L,$(_LIBDIRS))

_LIBS=m cisms
LIBS=$(addprefix -l,$(_LIBS))

CFLAGS=-O3 -Wall $(INCDIRS) -fopenmp
LDFLAGS=-O3 $(INCDIRS) $(LIBDIRS) $(LIBS) -Wl,-rpath,$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))/../lib -fopenmp

_SRCFILES=constrained_square.c
SRCS=$(_SRCFILES)

_OBJFILES=$(_SRCFILES:%.c=%.o)
OBJS=$(_OBJFILES)

EXES=constrained_square.exe

LIBNAME=libcisms.a

DIRGUARD=@mkdir -p $(@D)

all: $(EXES)

%.exe: $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

%.o : %.c
	$(CC) $< -c -o $@ $(CFLAGS)

.SECONDARY: $(OBJS)
.PHONY: clean

clean:
	rm *.o
	rm *.exe
