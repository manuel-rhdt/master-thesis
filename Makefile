
include Makefile.powerpc

DEST		= .

LINKER		= $(ARCHLINKER)

LDFLAGS 	= $(ARCHLDFLAGS) 

LIBS    	= $(ARCHLIBS)

DEBUG_CFLAGS	= $(ARCHDEBUGFLAGS) $(COMMON_CFLAGS) -DDEBUG

PROFILE_CFLAGS	= $(ARCHPROFFLAGS)

OPENGL_CFLAGS	= $(DEBUG_CFLAGS) -DUSE_OPENGL_VIEWER

COMMON_CFLAGS	= -Inumtools

CFLAGS		= $(ARCHCFLAGS) $(COMMON_CFLAGS)

CCFLAGS		= $(ARCHCCFLAGS)

CC		= $(ARCHCC)
CCC		= $(ARCHCCC)

SRCS		= main.c propagate.c read_comp_react.c\
		  stats.c analyse.c ran.n.c\
	 	  numtools/vectools.c numtools/numtools.c


OBJS		= $(ARCHOBJS) \
	          main.o propagate.o read_comp_react.o\
		  stats.o analyse.o ran.n.o\
	          numtools/vectools.o numtools/numtools.o

OUTPUT		= Gillespie.exe
OUTPUT2		= Gillespie.2.exe

vers2:  $(OBJS)
	$(LINKER) $(OBJS) $(LDFLAGS) $(LIBS) -o $(OUTPUT2)

all:	$(OBJS)
	$(LINKER) $(OBJS) $(LDFLAGS) $(LIBS) -o $(OUTPUT)

opengl:
	$(MAKE) "CFLAGS=$(OPENGL_CFLAGS)" all

debug:
	$(MAKE) "CFLAGS=$(DEBUG_CFLAGS)" all

profile:
	$(MAKE) "CFLAGS=-g -pg -Inumtools" pall

pall:	$(OBJS)
	$(LINKER) $(OBJS) $(LDFLAGS) $(LIBS) -g -pg -o $(OUTPUT)

depend: 
	makedepend -f Makefile.powerpc -- $(CFLAGS) -- $(SRCS)

clean:;	@echo "Cleaning... "
	@rm -f *.o numtools/*.o  "#"* core *.bak
	@echo "Done."



# DO NOT DELETE
