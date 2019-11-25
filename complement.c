/**
    This file is part of exactcolors.

    exactcolors is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    exactcolors is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with exactcolors.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <getopt.h>
#include <string.h>


#include "color.h"
#include "lp.h"
#include "graph.h"
#include "mwis.h"

#define STABLE_EPS  0.000001
#define STABLE_VIOL 0.001

#define STABLE_CLIQUE  1
#define STABLE_ODDHOLE 2

typedef struct branchobj {
    int vertex;
    int value;
} branchobj;

typedef struct stablecut {
    int count;
    int *inodes;
    int rhs;
    int type;
    struct stablecut *next;
} stablecut;

typedef struct prob {
    int *icuts;
    int cutcount;
    int cutspace;
    int branchdepth;
    branchobj *branchhistory;
    COLORlp_warmstart *warmstart;
    struct prob *next;
    struct prob *prev;
} prob;

typedef struct pool {
    stablecut *cuts;
    int cutcount;
    int cutspace;
} pool;

static char *graphfile = (char *) NULL;
static char *outfile = (char *) NULL;
static int debug = 0;
static int rundfs = 0;
static int usecuttingplanes = 0;
static int useostergard = 0;

int main (int ac, char **av);
static int parseargs (int ac, char **av);
static void usage (char *f);
static int get_problem_name(char* pname,const char* efname);

int main (int ac, char **av)
{
    char pname[256] = "";
    int rval = 0,i;
    int ncount, ecount;
    int *elist = (int *) NULL;
    int *wlen = (int *) NULL;
    COLORadjgraph G,Gc;

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    get_problem_name (pname, graphfile);

    rval = COLORread_dimacs (graphfile, &ncount, &ecount, &elist, &wlen);
    COLORcheck_rval (rval, "COLORread_dimacs failed");

    COLORadjgraph_init(&G);
    COLORadjgraph_init(&Gc);

    rval = COLORadjgraph_build(&G,ncount,ecount, elist);
    COLORcheck_rval (rval, "COLORadjgraph_build failed");

    rval =  COLORcheck_connectedness(&G);
    COLORcheck_rval (rval, "COLORcheck_connectedness failed");

    rval = COLORadjgraph_build_complement(&Gc,&G);
    COLORcheck_rval (rval, "COLORadjgraph_build_complement failed");

    rval =  COLORadjgraph_simplify(&Gc);
    COLORcheck_rval(rval,"COLORadjgraph_simplify");

    COLOR_IFFREE(elist,int);
    rval = COLORadjgraph_extract_edgelist(&ecount,&elist,&Gc);
    COLORcheck_rval (rval, "COLORadjgraph_extract_edgelist failed");

    printf("p edge %d %d\n",ncount, ecount);
    for (i = 0; i < ecount; ++i) {
       printf("e %d %d\n", elist[2*i]+1, elist[2*i+1] + 1);
    }

CLEANUP:
    COLORadjgraph_free(&G);
    COLORadjgraph_free(&Gc);

    COLOR_IFFREE(elist,int);
    COLOR_IFFREE(wlen,int);
    return rval;
}
static int parseargs (int ac, char **av)
{
    int c;
    int rval = 0;

    while ((c = getopt (ac, av, "bcdOo:")) != EOF) {
        switch (c) {
        case 'b':
            rundfs = 1;
            break;
        case 'c':
            usecuttingplanes = 1;
            break;
        case 'd':
            debug = 1;
            break;
        case 'O':
            useostergard = 1;
            break;
        case 'o':
            outfile = optarg;
            break;
        default:
            rval = 1;
            goto CLEANUP;
        }
    }

    if (ac <= optind) {
        rval = 1; goto CLEANUP;
    } else {
        graphfile = av[optind++];
    }

CLEANUP:

    if (rval) { usage (av[0]);}
    return rval;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage %s: [-see below-] graph_file\n", f);
    fprintf (stderr, "   -b    run DFS branching for stable sets\n");
    fprintf (stderr, "   -c    use stable-set cutting-plane algorithm\n");
    fprintf (stderr, "   -d    turn on debugging\n");
    fprintf (stderr, "   -O    Ostergard alg for max-weight clique\n");
    fprintf (stderr, "   -o f  write result file f\n");
    fprintf (stderr,"    NOTE: clique is found if -c not specified\n");
}

static int get_problem_name(char* pname,const char* efname)
{
    int rval = 0;
    int len = 0;
    const char *fname = strrchr(efname,'/');
    char *lastdot = strrchr(efname,'.');

    if(!fname) {
        /* no slashes in efname.*/
        fname = efname;
    } else {
        fname++;
    }

    if (lastdot) {
       len = lastdot - fname + 1;
    } else {
       len = strlen(fname);
    }

    if (snprintf(pname,len,"%s",fname) < 0) {
        fprintf (stderr, "snprintf failed\n");
        rval = 1; goto CLEANUP;
    }
    printf("c Extracted problem name %s\n",pname);  fflush (stdout);

CLEANUP:

   return rval;
}
