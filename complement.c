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
#include <fenv.h>

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

    rval = COLORadjgraph_build_complement(&Gc,&G);
    COLORcheck_rval (rval, "COLORadjgraph_build_complement failed");

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
            usage (av[0]);
            rval = 1;  goto CLEANUP;
        }
    }

    if (ac <= optind) {
        rval = 1; goto CLEANUP;
    } else {
        graphfile = av[optind++];
    }

CLEANUP:

    if (rval) usage (av[0]);
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
    printf("Extracted problem name %s\n",pname);  fflush (stdout);

CLEANUP:

   return 0;
}

#if 0
static int read_dimacs_mwis (char *f, int *pncount, int *pecount, int **pelist,
        int **pwlen)
{
    int rval = 0;
    int ncount, ecount, icount = 0, haveprob = 0;
    int end0, end1, n, i, len;
    double t;
    int *elist = (int *) NULL;
    int *wlen = (int *) NULL;
    char buf[256], *p;
    FILE *in = (FILE *) NULL;

    in = fopen (f, "r");
    if (!in) {
        fprintf (stderr, "Unable to open %s for input\n", f);
        rval = 1;  goto CLEANUP;
    }

    while (fgets (buf, 254, in) != (char *) NULL) {
        p = buf;
        if (p[0] == 'c') {
            printf ("Comment: %s", p+1);
        } else if (p[0] == 'p') {
            const char* delim = " \t\n";
            char* data = (char *) NULL;
            if (haveprob) {
                fprintf (stderr, "ERROR in Dimacs file -- two p lines\n");
                rval = 1;  goto CLEANUP;
            }
            haveprob = 1;
            data = strtok(p,delim); /* get 'p' */

            data = strtok(NULL,delim); /* get type */
            if ( strcmp(data,"graph") && strcmp(data,"edges") &&
                                        strcmp(data,"col") ) {
                fprintf (stderr, "ERROR in Dimacs file -- not a graph file\n");
                rval = 1;  goto CLEANUP;
            }
            data = strtok(NULL,delim);
            sscanf (data, "%d", &ncount);
            data = strtok(NULL,delim);
            sscanf (data, "%d", &ecount);

            printf ("Number of Nodes: %d\n", ncount);
            printf ("Number of Edges: %d\n", ecount);
            elist = COLOR_SAFE_MALLOC (2 * ecount, int);
            COLORcheck_NULL(elist, "out of memory for elist");
            wlen = COLOR_SAFE_MALLOC (ncount, int);
            COLORcheck_NULL(elist, "out of memory for wlen");
            for (i = 0; i < ncount; i++) wlen[i] = 0;
        } else if (p[0] == 'e') {
            if (!haveprob) {
                fprintf (stderr, "ERROR in Dimacs file -- e before p\n");
                rval = 1;  goto CLEANUP;
            }
            if (icount >= ecount) {
                fprintf (stderr, "ERROR in Dimacs file -- to many edges\n");
                rval = 1;  goto CLEANUP;
            }
            p++;
            sscanf (p, "%d %d", &end0, &end1);
            elist[2*icount] = end0-1;    /* Number nodes from 0, not 1 */
            elist[2*icount+1] = end1-1;
            icount++;
        } else if (p[0] == 'n') {
            if (!haveprob) {
                fprintf (stderr, "ERROR in Dimacs file -- n before p\n");
                rval = 1;  goto CLEANUP;
            }
            p++;
            sscanf (p, "%d %d", &n, &len);
            wlen[n-1] = len;
        }
    }

    *pncount = ncount;
    *pecount = icount;
    *pelist = elist;
    *pwlen = wlen;

CLEANUP:

    if (rval) {
        COLOR_IFFREE(elist,int);
        COLOR_IFFREE(wlen,int);
    }
    if (in) fclose (in);

    return rval;
}
#endif

static int build_set (int *nlist, int count, COLORset **pcset)
{
    int i, rval = 0;
    COLORset *cset = (COLORset *) NULL;

    cset = COLOR_SAFE_MALLOC (1, COLORset);
    COLORcheck_NULL (cset, "out of memory for cset");
    COLORinit_set (cset);

    cset->members = COLOR_SAFE_MALLOC (count, int);
    COLORcheck_NULL (cset->members, "out of memory for members");
    for (i = 0; i < count; i++) cset->members[i] = nlist[i];
    cset->count = count;

    *pcset = cset;

CLEANUP:

    if (rval) {
        COLORfree_set (cset);
        COLOR_IFFREE (cset, COLORset);
    }

    return rval;
}

static int build_stablecut (int *nlist, int count, int cuttype,
        stablecut **pcut)
{
    int i, rval = 0;
    stablecut *cut = (stablecut *) NULL;

    cut = COLOR_SAFE_MALLOC (1, stablecut);
    COLORcheck_NULL (cut, "out of memory for cut");
    init_stablecut (cut);

    cut->inodes = COLOR_SAFE_MALLOC (count, int);
    COLORcheck_NULL (cut->inodes, "out of memory for inodes");
    for (i = 0; i < count; i++) cut->inodes[i] = nlist[i];
    cut->count = count;
    cut->type = cuttype;
    if (cuttype == STABLE_CLIQUE) {
        cut->rhs = 1;
    } else if (cuttype == STABLE_ODDHOLE) {
        if (count % 2 == 0) {
            printf ("Odd hole with an even number of nodes\n");
            rval = 1; goto CLEANUP;
        }
        cut->rhs = count/2;
    } else {
        printf ("Unknown cut type: %d\n", cut->type);
        rval = 1;  goto CLEANUP;
    }

    *pcut = cut;

CLEANUP:

    if (rval) {
        free_stablecut (cut);
        COLOR_IFFREE (cut, stablecut);
    }

    return rval;
}

void init_stablecut (stablecut *cut)
{
    if (cut) {
        cut->count = 0;
        cut->inodes = (int *) NULL;
        cut->rhs = 0;
        cut->type = 0;
        cut->next = (stablecut *) NULL;
    }
}

void free_stablecut (stablecut *cut)
{
    if (cut) {
        COLOR_IFFREE (cut->inodes, int);
        init_stablecut (cut);
    }
}

static int call_lp_solver (COLORlp *lp, double *objval, double *x)
{
    int rval = 0;

    rval = COLORlp_optimize (lp);
    COLORcheck_rval (rval, "COLORlp_optimize failed");

    if (objval) {
        rval = COLORlp_objval (lp, objval);
        COLORcheck_rval (rval, "COLORlp_objval failed");
    }

    if (x) {
        rval = COLORlp_x (lp, x);
        COLORcheck_rval (rval, "COLORlp_x failed");
    }

CLEANUP:

    return rval;
}

static void perm_dbl_quicksort (int *perm, const double *len, int n)
{
    int i, j, temp;
    double t;

    if (n <= 1) return;

    COLOR_SWAP (perm[0], perm[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do i++; while (i < n && len[perm[i]] < t);
        do j--; while (len[perm[j]] > t);
        if (j < i) break;
        COLOR_SWAP (perm[i], perm[j], temp);
    }
    COLOR_SWAP (perm[0], perm[j], temp);

    perm_dbl_quicksort (perm, len, j);
    perm_dbl_quicksort (perm + i, len, n - i);
}

static void perm_dbl_rquicksort (int *perm, const double *len, int n)
{
    int i, j, temp;
    double t;

    if (n <= 1) return;

    COLOR_SWAP (perm[0], perm[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do i++; while (i < n && len[perm[i]] > t);
        do j--; while (len[perm[j]] < t);
        if (j < i) break;
        COLOR_SWAP (perm[i], perm[j], temp);
    }
    COLOR_SWAP (perm[0], perm[j], temp);

    perm_dbl_rquicksort (perm, len, j);
    perm_dbl_rquicksort (perm + i, len, n - i);
}

int COLORdbg_lvl() {
   return debug;
}

