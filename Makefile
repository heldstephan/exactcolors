# Adapt GUPATH to point to your gurobi installation
# or set the environment variable GUROBI_HOME accordingly
GUPATH=$(GUROBI_HOME)
QSPATH=/home/fac/bico/QS/work

#LPINCLUDE=$(GUPATH)/include
#LPLIB=$(GUPATH)/lib/libgurobi.so.2.0.1
#LPSOURCE=lpgurobi.o
#GRBMWIS=mwis_grb.o
#GUROBI_FLAG=-DUSE_GUROBI

LPINCLUDE=$(QSPATH)
LPLIB=$(QSPATH)/qsopt.a
LPSOURCE=lpqsopt.o


# SEWELL_FLAG=-DHAVE_SEWELL
# SEWELL_LIB=-L . -lsewell

CC=gcc
CFLAGS=  -O3  -std=c99 -pedantic -Wall -Wshadow -W -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wpointer-arith -Wnested-externs -Wundef -Wcast-qual -Wcast-align -Wwrite-strings -I$(LPINCLUDE) $(SEWELL_FLAG)
OBJFILES=color.o graph.o greedy.o $(LPSOURCE) mwis.o $(GRBMWIS) mwis_grdy.o plotting.o heap.o util.o cliq_enum.o bbsafe.o
STABFILES=stable.o graph.o greedy.o util.o $(LPSOURCE) cliq_enum.o
BOSSFILES=graph.o bbsafe.o util.o
CBOSSFILES=color_main.o $(OBJFILES)
CWORKERFILES=color_worker.o $(OBJFILES)

all: color color_worker stable queen test_boss test_worker test_tell

color: $(CBOSSFILES)
	$(CC) $(CFLAGS) -o color $(CBOSSFILES) $(LPLIB) -lm -lpthread $(SEWELL_LIB)

color_worker: $(CWORKERFILES)
	$(CC) $(CFLAGS) -o color_worker $(CWORKERFILES) $(LPLIB) -lm -lpthread $(SEWELL_LIB)

stable: $(STABFILES)
	$(CC) $(CFLAGS) -o stable $(STABFILES) $(LPLIB) -lm -lpthread

queen: queen.c
	$(CC) $(CFLAGS) -o queen queen.c -lm -lpthread

test_boss: test_boss.o $(BOSSFILES)
	$(CC) $(CFLAGS) -o test_boss test_boss.o $(BOSSFILES) -lm -lpthread

test_worker: test_worker.o $(BOSSFILES)
	$(CC) $(CFLAGS) -o test_worker test_worker.o $(BOSSFILES) -lm -lpthread

test_tell: test_tell.o $(BOSSFILES)
	$(CC) $(CFLAGS) -o test_tell test_tell.o $(BOSSFILES) -lm -lpthread

tags:
	etags *.[hc]
clean:
	rm -f *.o color stable test_boss test_worker test_tell mwis_gurobi.log gurobi.log look.lp vg.log*

color.o:     color_main.c color.c color.h color_private.h lp.h color_defs.h mwis.h plotting.h heap.h bbsafe.h
color_worker.o: color_worker.c color_private.h color_defs.h bbsafe.h
heap.o:      heap.c heap.h color_defs.h
graph.o:     graph.c graph.h color_defs.h
greedy.o:    greedy.c  color.h graph.h color_defs.h
lpgurobi.o:  lpgurobi.c color.h lp.h color_defs.h
lpqsopt.o:   lpqsopt.c color.h lp.h color_defs.h
mwis.o:      mwis.c mwis.h color.h color_defs.h
mwis_grdy.o: mwis_grdy.c color.h graph.h color_defs.h heap.h
mwis_grb.o:  mwis_grb.c color.h lp.h color_defs.h
stable.o:    stable.c color.h graph.h lp.h
util.o:      util.c color.h
cliq_enum.o: color.h lp.h graph.h mwis.h
test_boss.o: test_boss.c bbsafe.h
test_worker.o: test_worker.c bbsafe.h
test_tell.o: test_tell.c bbsafe.h

