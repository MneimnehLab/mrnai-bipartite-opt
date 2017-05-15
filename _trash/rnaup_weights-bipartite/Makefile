CC=gcc
CFLAGS=-c -Wall -w 
LDFLAGS=

WEIGHTS_RNAUP_SOURCE = weights_rnaup.c

RNAUP_SOURCES = rnauplib/alifold.c rnauplib/cofold.c rnauplib/dist_vars.c rnauplib/duplex.c rnauplib/energy_par.c  rnauplib/fold.c rnauplib/fold_vars.c  rnauplib/list.c  rnauplib/params.c rnauplib/part_func.c rnauplib/part_func_co.c rnauplib/part_func_up.c  rnauplib/ribo.c rnauplib/stringdist.c rnauplib/subopt.c rnauplib/utils.c rnauplib/RNAup.c 
RNAUP_OBJECTS = $(RNAUP_SOURCES:.c=.o)

BY_BOTH_SOURCES = extractor/datastruct.c extractor/parse.c extractor/precomp.c 
BY_BOTH_OBJECTS = $(BY_BOTH_SOURCES:.c=.o)

WEIGHTS_RNAUP_OBJECTS = $(WEIGHTS_RNAUP_SOURCE:.c=.o)

WEIGHTS_RNAUP_EXEC = weights_rnaup

all: weights_rnaup

weights_rnaup: $(WEIGHTS_RNAUP_EXEC)


$(WEIGHTS_RNAUP_EXEC): $(RNAUP_OBJECTS) $(BY_BOTH_OBJECTS) $(WEIGHTS_RNAUP_OBJECTS)
	$(CC) $(LDFLAGS)   $(RNAUP_OBJECTS) $(BY_BOTH_OBJECTS) $(WEIGHTS_RNAUP_OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm -f *.o
	rm -f rnauplib/*.o
	rm -f extractor/*.o