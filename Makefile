CC=gcc

_INCDIRS=include
INCDIRS=$(addprefix -I,$(_INCDIRS))

_LIBDIRS=lib
LIBDIRS=$(addprefix -L,$(_LIBDIRS))

_LIBS=m
LIBS=$(addprefix -l,$(_LIBS))

CFLAGS=-O3 -Wall $(INCDIRS)
LDFLAGS=-O3 -shared $(INCDIRS) $(LIBDIRS) $(LIBS)

SRCDIR=src
_SRCFILES=particle.c config.c
SRCS=$(addprefix $(SRCDIR)/,$(_SRCFILES))

OBJDIR=build
_OBJFILES=$(_SRCFILES:%.c=%.o)
OBJS=$(addprefix $(OBJDIR)/,$(_OBJFILES))

LIBDIR=lib
LIBNAME=libcisms.so

ROOTDIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
INSTALLDIR=/usr/local/lib
INSTALLNAME=libcisms.so

DIRGUARD=@mkdir -p $(@D)

all: $(LIBDIR)/$(LIBNAME)

$(LIBDIR)/$(LIBNAME): $(OBJS)
	$(DIRGUARD)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)
	@ln -s $(ROOTDIR)/$(LIBDIR)/$(LIBNAME) $(INSTALLDIR)/$(INSTALLNAME)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(DIRGUARD)
	$(CC) $< -c -o $@ $(CFLAGS)

install:
	@ln -s $(LIBDIR)/$(LIBNAME) $(INSTALLDIR)/$(INSTALLNAME)

.SECONDARY: $(OBJS)
.PHONY: clean

clean:
	rm $(OBJDIR)/*.o
	rm $(LIBDIR)/$(LIBNAME)
	rm $(INSTALLDIR)/$(INSTALLNAME)
	rmdir $(OBJDIR) $(LIBDIR)
