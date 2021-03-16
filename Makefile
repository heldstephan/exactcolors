#     This file is part of exactcolors.

#     exactcolors is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     exactcolors is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with exactcolors.  If not, see <http://www.gnu.org/licenses/>.


# exactcolors requires one of the following 3 LP solvers: Gurobi, Cplex, or QSopt.
# Please set the environment path and uncomment the lines correspondingly.
# You might also need to adopt the  LPINCLUDE & LPLIB  paths further below.
GUPATH=$(GUROBI_HOME)
CPLEXPATH=$(CPLEX_HOME)
QSPATH=$(QSOPT_HOME)

USE_UBSAN=0
#CFLAGS+= -g
CFLAGS+= -O3


ifneq ($(QSPATH),)
LPINCLUDE=$(QSPATH)
LPLIB=$(QSPATH)/qsopt.a
LPSOURCE=lpqsopt.o
endif

ifneq ($(GUPATH),)
LPINCLUDE=$(GUPATH)/include
LPLIB=$(GUPATH)/lib/libgurobi*.so
LPSOURCE=lpgurobi.o
GRBMWIS=mwis_grb.o
GUROBI_FLAG=-DUSE_GUROBI
endif

ifneq ($(CPLEXPATH),)
PROCESSOR := $(shell uname -p)
LPINCLUDE=$(CPLEXPATH)/include/ilcplex
LPLIB=$(CPLEXPATH)/lib/x86-64_linux/static_pic/libcplex.a
ifeq ($(PROCESSOR), i686)
LPLIB=$(CPLEXPATH)/lib/x86_linux/static_pic/libcplex.a
endif
LPSOURCE=lpcplex.o
GRBMWIS=
GUROBI_FLAG=
endif


export CC=gcc
# Alternative compiler, e.g. for static code analysis.
#export CC=clang

export LD=gcc



#
# Valgrind does not support fegetround & fesetround. With following compile option
# their use is circumvented. We also recommend to use QSopt as the LP-solver while
# debugging with valgrind, as the commercial solvers impose valgrind errors internally.
#
#CFLAGS+= -DCOMPILE_FOR_VALGRIND



####################################################
# Below this comment changes should be unnecessary.#
####################################################

SEWELL_DIR=mwis_sewell
SEWELL_LDFLAG=-L $(SEWELL_DIR) -lsewell
SEWELL_LIB=$(SEWELL_DIR)/libsewell.a

EXACTCOLOR_DIR=.
EXACTCOLOR_LDFLAG=-L$(EXACTCOLOR_DIR) -lexactcolor $(LPLIB)  -ldl -lm -lpthread
EXACTCOLOR_LIB= $(EXACTCOLOR_DIR)/libexactcolor.a 

CFLAGS += -std=c99 -D_XOPEN_SOURCE=700 -pedantic -Wall -Wshadow -W -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wpointer-arith -Wnested-externs -Wundef -Wcast-qual -Wcast-align -Wwrite-strings -I$(LPINCLUDE)
export CFLAGS

# UBSAN
ifeq ($(USE_UBSAN), 1)
	CFLAGS += -fsanitize=undefined -fsanitize=float-divide-by-zero
endif


OBJFILES=color.o color_backup.o color_parms.o graph.o greedy.o $(LPSOURCE) mwis.o $(GRBMWIS) mwis_grdy.o plotting.o heap.o util.o cliq_enum.o bbsafe.o rounding_mode.o
STABFILES=stable.o graph.o greedy.o util.o $(LPSOURCE) cliq_enum.o
STABGRDYFILES=stable_grdy.o graph.o greedy.o util.o $(LPSOURCE) cliq_enum.o mwis.o mwis_grdy.o  heap.o $(SEWELL_LIB)
BOSSFILES=graph.o bbsafe.o util.o rounding_mode.o
CBOSSFILES=color_version.h color_main.o rounding_mode.o
CWORKERFILES=color_worker.o
CKILLERFILES=color_jobkiller.o
PARTFILES=partition.o
COMPFILES=complement.o

all: color color_worker color_jobkiller stable stable_grdy queen test_boss test_worker test_tell partition complement dsatur


scan_build: *.[hc] mwis_sewell/*.[hc]
	export CC=ccc-analyzer
	scan-build -v -o clang make -j


testmyciel4:
	./color test/instances/myciel4.col  |grep LB > test/myciel4.con
	diff test/myciel4.con test/golden/myciel4.color.con
	./dsatur test/instances/myciel4.col  |grep LB > test/myciel4.con
	diff test/myciel4.con test/golden/myciel4.dsatur.con

testqueen8:
	./color test/instances/queen8_8.col  |grep LB > test/queen8_8.con
	diff test/queen8_8.con test/golden/queen8_8.color.con
	./dsatur test/instances/queen8_8.col  |grep LB > test/queen8_8.con
	diff test/queen8_8.con test/golden/queen8_8.dsatur.con

test: testmyciel4 testqueen8



libexactcolor.a: $(OBJFILES)
	$(AR) rcs libexactcolor.a $(OBJFILES)

color: $(EXACTCOLOR_LIB) $(SEWELL_LIB) $(CBOSSFILES) color_worker
	$(LD)  $(CFLAGS)  -o color color_main.o $(EXACTCOLOR_LDFLAG) $(SEWELL_LDFLAG)  

color_worker: $(EXACTCOLOR_LIB) $(SEWELL_LIB) $(CWORKERFILES)
	$(CC) $(CFLAGS) -o color_worker $(CWORKERFILES)  $(EXACTCOLOR_LDFLAG) $(SEWELL_LDFLAG) 

color_jobkiller:  $(EXACTCOLOR_LIB) $(SEWELL_LIB) $(CKILLERFILES)
	$(CC) $(CFLAGS) -o color_jobkiller $(CKILLERFILES)  $(EXACTCOLOR_LDFLAG) $(SEWELL_LDFLAG)

$(SEWELL_LIB): $(SEWELL_DIR)/*[hc] $(SEWELL_DIR)/Makefile
	cd $(SEWELL_DIR) && $(MAKE) USE_UBSAN=$(USE_UBSAN)

stable: $(EXACTCOLOR_LIB) $(STABFILES)
	$(CC) $(CFLAGS) -o stable $(STABFILES)  $(EXACTCOLOR_LDFLAG)

stable_grdy: $(EXACTCOLOR_LIB) $(STABGRDYFILES)
	$(CC) $(CFLAGS) -o stable_grdy $(STABGRDYFILES)  $(EXACTCOLOR_LDFLAG)

partition: $(EXACTCOLOR_LIB) $(SEWELL_LIB) $(PARTFILES)
	$(CC) $(CFLAGS) -o partition $(PARTFILES) $(EXACTCOLOR_LDFLAG) $(SEWELL_LIB)

complement: $(EXACTCOLOR_LIB) $(SEWELL_LIB) $(COMPFILES)
	$(CC) $(CFLAGS) -o complement $(COMPFILES)  $(EXACTCOLOR_LDFLAG) $(SEWELL_LIB)

dsatur: dsatur.o graph.o color.o rounding_mode.o $(EXACTCOLOR_LIB) $(SEWELL_LIB)
	$(LD) $(CFLAGS) -o dsatur dsatur.o graph.o color.o color_parms.o rounding_mode.o $(EXACTCOLOR_LDFLAG) $(SEWELL_LDFLAG)

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
	rm -f *.o color stable test_boss test_worker test_tell partition mwis_gurobi.log gurobi.log look.lp vg.log* color_version.h color_worker color_jobkiller queen complement stable_grdy libexactcolor.a
	rm -rf clang
	rm -f libexactcolor
	rm -f test/*.con
	cd $(SEWELL_DIR) && $(MAKE) clean

SRCFILES=bbsafe.c color_backup.c  color.h color_parms.c  graph.c   heap.c     lpgurobi.c  mwis.c       mwis.h                  mwis_sewell/mwss_ext.h  partition.c  queen.c        test_boss.c    util.c bbsafe.h     color.c         color_jobkiller.c  color_parms.h    color_worker.c   graph.h   heap.h     lp.h        mwis_grb.c   mwis_sewell/mwss.c      mwis_sewell/mwss.h      plotting.c   stable.c       test_tell.c cliq_enum.c  color_defs.h    color_main.c       color_private.h  complement.c     greedy.c  lpcplex.c  lpqsopt.c   mwis_grdy.c  mwis_sewell/mwss_ext.c  mwis_sewell/wstable.c   plotting.h   stable_grdy.c  test_worker.c

color_version.h: $(SRCFILES)
	./create_version_header > color_version.h

color.o:     color_main.c color.c color.h color_private.h lp.h color_defs.h mwis.h plotting.h heap.h bbsafe.h color_version.h
color_worker.o: color_worker.c color_private.h color_defs.h bbsafe.h
color_backup.o: color_backup.c color_private.h color_defs.h
color_parms.o: color_parms.c color_parms.h color_defs.h
partition.o: partition.c  color.h graph.h color_defs.h
complement.c:  color.h graph.h color_defs.h
dsatur.o:    dsatur.c graph.h color_defs.h color.h color_private.h color_parms.h rounding_mode.h
heap.o:      heap.c heap.h color_defs.h
graph.o:     graph.c graph.h color_defs.h
greedy.o:    greedy.c  color.h graph.h color_defs.h
lpgurobi.o:  lpgurobi.c color.h lp.h color_defs.h
lpcplex.o:   lpcplex.c color.h lp.h color_defs.h
lpqsopt.o:   lpqsopt.c color.h lp.h color_defs.h
mwis.o:      mwis.c mwis.h color.h color_defs.h
mwis_grdy.o: mwis_grdy.c color.h graph.h color_defs.h heap.h
mwis_grb.o:  mwis_grb.c color.h lp.h color_defs.h
stable.o:    stable.c color.h graph.h lp.h
util.o:      util.c color.h
cliq_enum.o: color.h lp.h graph.h mwis.h
rounding_mode.o: rounding_mode.c rounding_mode.h
test_boss.o: test_boss.c bbsafe.h
test_worker.o: test_worker.c bbsafe.h
test_tell.o: test_tell.c bbsafe.h
