include ../Makefile.in

BINDIR = ..

INCLUDES = -I/p/include -I../multilevelLib -I../metisLib
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)

LIBS = -lmultilevel -lmetis -lm

LD = $(CC) $(LDOPTIONS) -L/p/lib -L. -L..
METISLIB = ../libmetis.a ../libmultilevel.a

MLKKMOBJS = neo.o io.o


.c.o:
	$(CC) $(CFLAGS) -c $*.c

.cpp.o:
	$(CC) $(CFLAGS) -c $*.cpp

default: $(BINDIR)/neo$(VERNUM)


$(BINDIR)/neo$(VERNUM): $(MLKKMOBJS) $(METISLIB)
	$(LD)  -o $@ $(CFLAGS) $(MLKKMOBJS) $(LIBS)
	chmod 755 $@



clean:
	rm -f *.o

realclean:
	rm -f *.o ;\
	rm -f $(BINDIR)/neo$(VERNUM) 
