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
static int optimal_stable_set (int ncount, int ecount, int *elist,
    int *weights);
static int cutting_loop (COLORadjgraph *G, COLORlp *lp, pool *P,
    prob *probdata, int silent);
static int add_cuts (COLORlp *lp, pool *P, prob *sp, stablecut **pclist);
static int build_stable_lp (COLORlp **lp, COLORadjgraph *G, int *weights,
    int *clist);
static int add_cut_to_lp (COLORlp *lp, stablecut *cset);
static void init_stable_prob (prob *p);
static void free_stable_prob (prob *p);
static int build_stable_prob (prob *root, int setcount, int ncount);
static int add_cut_to_prob (prob *p, int ind);
static void init_cut_pool (pool *P);
static int build_cut_pool (pool *P, int setcount, int ncount);
static int add_cut_to_pool (pool *P, stablecut **pcset, int *hit);
static void free_cut_pool (pool *P);
static int cover_edge_cliques (COLORadjgraph *G, int ecount, int *elist,
    stablecut **pclist, int *setcount);
static int find_violated_cliques (COLORadjgraph *G, double *x, 
    stablecut **pclist, int *pcount);
static int cover_node (COLORadjgraph *G, double *weights, int k, 
    stablecut **pcset, int *nmarks, int *nlist);
static int cover_edge (COLORadjgraph *G, int end1, int end2, stablecut **pcset,
    int *nmarks, int *nlist);
static void mark_edges (COLORadjgraph *G, int **incid, stablecut *cset, 
    int *emarks, int *nmarks, int *elist);
static void check_clique (COLORadjgraph *G, int *nmarks, stablecut *cset, 
    int *yesno);
static int find_violated_holes (COLORadjgraph *G, double *x,
    stablecut **cutlist, int *pcount);
static int stable_greedy_x (COLORadjgraph *G, int *weights, double *x,
    COLORset **pcset, int *wt);
static int call_dfs_branching (COLORadjgraph *G, int *weights, COLORlp *lp,
    pool *P, prob *probdata, int *bestval);
static int dfs_branching (COLORadjgraph *G, int *weights, COLORlp *lp, pool *P,
    prob *probdata, int *bestval, double *x, int depth, int *bcount);
static int find_branch (COLORadjgraph *G, double *x, int *pvertex);
static int get_integer_soln (COLORadjgraph *G, int *weights, double *x, 
    COLORset *sol, int *solval);
static int build_incidence (COLORadjgraph *G, int ecount, int *elist, 
    int **incid);
static int parseargs (int ac, char **av);
static void usage (char *f);
static int get_problem_name(char* pname,const char* efname);
static void perm_dbl_quicksort (int *perm, const double *len, int n);
static void perm_dbl_rquicksort (int *perm, const double *len, int n);
static int build_set (int *nlist, int count, COLORset **pcset);
void init_stablecut (stablecut *cut);
void free_stablecut (stablecut *cut);
static int build_stablecut (int *nlist, int count, int cuttype,
        stablecut **pcut);
static int call_lp_solver (COLORlp *lp, double *objval, double *x);

int main (int ac, char **av)
{
    char pname[256] = "";
    int rval = 0;
    int ncount, ecount, val, i;
    int *elist = (int *) NULL;
    int *wlen = (int *) NULL;
    COLORset *cliques = (COLORset *) NULL;
    int ncliques = 0;
    double szeit;
    int    nrbranches = INT_MAX;
    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    get_problem_name (pname, graphfile);

    rval = COLORread_dimacs (graphfile, &ncount, &ecount, &elist, &wlen);
    COLORcheck_rval (rval, "COLORread_dimacs failed");

    szeit = COLORutil_zeit ();

    if (usecuttingplanes) {
        rval = optimal_stable_set (ncount, ecount, elist, wlen);
        COLORcheck_rval (rval, "optimal_stable_set failed");
    } else {
        if (useostergard) {
            rval = COLORclique_ostergard (&cliques, &ncliques, ncount, ecount,
                                          elist, wlen, COLOR_MAXINT, &val,
                                          nrbranches);
            COLORcheck_rval (rval, "COLORcliq_ostergard failed");
        } else {
            rval = COLORclique_enum (&cliques, &ncliques, ncount, ecount, elist,
                                     wlen, COLOR_MAXINT, &val);
            COLORcheck_rval (rval, "COLORcliq_enum failed");
        }
        
        printf ("Optimal Weight Clique: %d\n", val);
        fflush (stdout);

        printf ("Clique: ");
        for (i = 0; i < cliques[0].count; i++) {
            printf ("%d ", cliques[0].members[i]);
       }
       printf ("\n");
       COLORfree_set (cliques);
       COLOR_FREE (cliques, COLORset);
    }

    printf ("Running Time: %.2lf seconds\n", COLORutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    COLOR_IFFREE(elist,int);
    COLOR_IFFREE(wlen,int);
    return rval;
}

static int optimal_stable_set (int ncount, int ecount, int *elist, int *weights)
{
    int rval = 0;
    stablecut *clist = (stablecut *) NULL;
    COLORset *gcset = (COLORset *) NULL;
    int gwt = 0, setcount = 0, bestweight = 0;
    COLORadjgraph G;
    pool P;
    prob *root = (prob *) NULL;
    COLORlp *lp = (COLORlp *) NULL;
    double val, gap;
    double *x = (double *) NULL;

    init_cut_pool (&P);
    rval = COLORadjgraph_build (&G, ncount, ecount, elist);
    COLORcheck_rval (rval, "COLORadjgraph_build failed");

    root = COLOR_SAFE_MALLOC (1, prob);
    COLORcheck_NULL (root, "out of memory for root");
    init_stable_prob (root);

    rval = build_stable_lp (&lp, &G, weights, (int *) NULL);
    COLORcheck_rval (rval, "build_stable_lp failed");

    rval = cover_edge_cliques (&G, ecount, elist, &clist, &setcount);
    COLORcheck_rval (rval, "cover_edge_cliques");

    rval = build_cut_pool (&P, setcount, G.ncount); 
    COLORcheck_rval (rval, "build_cut_pool failed");

    rval = build_stable_prob (root, setcount, G.ncount);
    COLORcheck_rval (rval, "build_stable_prob failed");

    rval = add_cuts (lp, &P, root, &clist);
    COLORcheck_rval (rval, "add_cuts failed");

    x = COLOR_SAFE_MALLOC (ncount, double);
    COLORcheck_NULL (x, "out of memory for x");

    rval = call_lp_solver (lp, &val, (double *) NULL);
    COLORcheck_rval (rval, "call_lp_solver failed");

    printf ("Initial LP Objective: %.6lf\n", val);
    fflush (stdout); 

    rval = cutting_loop (&G, lp, &P, root, 0);
    COLORcheck_rval (rval, "cutting_loop failed");
    printf ("Total Cuts: %d\n", root->cutcount); 

    rval = call_lp_solver (lp, &val, x);
    COLORcheck_rval (rval, "call_lp_solver failed");

    printf ("Call stable_greedy_x ...\n"); fflush (stdout);

    rval = stable_greedy_x (&G, weights, x, &gcset, &gwt);
    COLORcheck_rval (rval, "stable_greedy_x failed");
    printf ("Initial Stable Set: %d\n", gwt); fflush (stdout);

    gap = val - ((double) gwt);
    if (gap < 1.0-STABLE_EPS) {
        printf ("Optimized at Root.\n");
        goto CLEANUP;
    }

    printf ("Root Optimality Gap: %.6lf\n", gap);
    fflush (stdout);

    if (rundfs) {
        bestweight = gwt;
        rval = call_dfs_branching (&G, weights, lp, &P, root, &bestweight);
        COLORcheck_rval (rval, "call_dfs_branching failed");
    }

CLEANUP:

    if (lp) COLORlp_free (&lp);
    if (root) {
        free_stable_prob (root);
        COLOR_FREE (root, prob);
    }
    if (gcset) {
        COLORfree_set (gcset);
        COLOR_FREE (gcset, COLORset);
    }
    COLOR_IFFREE (x, double);
    COLORadjgraph_free (&G);
    free_cut_pool (&P);
    return rval;
}

static int cutting_loop (COLORadjgraph *G, COLORlp *lp, pool *P,
        prob *probdata, int silent)
{
    int rval = 0;
    int setcount = 0, iterations = 0, hcount = 0;
    double val;
    double *x = (double *) NULL;
    stablecut *clist, *hlist;

    x = COLOR_SAFE_MALLOC (G->ncount, double);
    COLORcheck_NULL (x, "out of memory for x");

    rval = call_lp_solver (lp, &val, x);
    COLORcheck_rval (rval, "call_lp_solver failed");

    iterations = 0;
    do {
        setcount = 0; hcount = 0;
        clist = (stablecut *) NULL;
        rval = find_violated_cliques (G, x, &clist, &setcount);
        COLORcheck_rval (rval, "find_violated_cliques failed");

        if (setcount > 0) {
            rval = add_cuts (lp, P, probdata, &clist);
            COLORcheck_rval (rval, "add_cuts failed");

            rval = call_lp_solver (lp, &val, x);
            COLORcheck_rval (rval, "call_lp_solver failed");

            if (!silent) {
                printf ("Iteration %d: added %2d cliques,   LP = %.6lf\n",
                         iterations, setcount, val);
                fflush (stdout);
            }
        }

        if (setcount <= 1) {
            hlist = (stablecut *) NULL;
            rval = find_violated_holes (G, x, &hlist, &hcount);
            COLORcheck_rval (rval, "find_violated_holes failed");

            if (hcount > 0) {
                rval = add_cuts (lp, P, probdata, &hlist);
                COLORcheck_rval (rval, "add_cuts failed");

                rval = call_lp_solver (lp, &val, x);
                COLORcheck_rval (rval, "call_lp_solver failed");

                if (!silent) {
                    printf ("Iteration %d: added %2d odd holes, LP = %.6lf\n",
                             iterations, hcount, val);
                    fflush (stdout);
                }
            }
        }
        iterations++;
    } while ((setcount > 0 || hcount > 0) && iterations < 500);

CLEANUP: 

    COLOR_IFFREE (x, double);
    return rval;
}

static int add_cuts (COLORlp *lp, pool *P, prob *sp, stablecut **pclist)
{
    int hit, rval = 0;
    stablecut *clist = *pclist;
    stablecut *cnext;
    
    while (clist) {
        cnext = clist->next;
        if (clist->type != STABLE_CLIQUE && clist->type != STABLE_ODDHOLE) {
            printf ("Unknown cut type: %d\n", clist->type);
            rval = 1; goto CLEANUP;
        }
        rval = add_cut_to_lp (lp, clist);
        COLORcheck_rval (rval, "add_to_lp failed");
        rval = add_cut_to_pool (P, &clist, &hit);
        COLORcheck_rval (rval, "add_cut_to_pool failed");
        rval = add_cut_to_prob (sp, hit);
        COLORcheck_rval (rval, "add_cut_to_prob failed");
        clist = cnext;
    }
    *pclist = (stablecut *) NULL;

CLEANUP:

    return rval;
}

static int build_stable_lp (COLORlp **lp, COLORadjgraph *G, int *weights,
        int *clist)
{
    int rval = 0;
    int i, ncount = G->ncount;
    double w;

    rval = COLORlp_init (lp, "stableset");
    COLORcheck_rval (rval, "COLORlp_init failed");
    rval = COLORlp_objective_sense (*lp, COLORlp_MAX);  

    for (i = 0; i < ncount; i++) {
        w = (double) weights[i];
        rval = COLORlp_addcol (*lp, 0, (int *) NULL, (double *) NULL,
                 w, 0.0, 1.0, COLORlp_CONTINUOUS, NULL);
        COLORcheck_rval (rval, "COLORlp_addcol failed");
    }

    if (clist) {
        printf ("ERROR: Not yet set up for clique list\n");
        rval = 1;  goto CLEANUP;
    }

CLEANUP:

    return rval;
}

static int add_cut_to_lp (COLORlp *lp, stablecut *c) 
{
    int rval = 0;
    int i;
    double rhs;
    double *coef = (double *) NULL;

    if (c->type != STABLE_CLIQUE && c->type != STABLE_ODDHOLE) {
        printf ("Unknon cut type: %d\n", c->type);
        rval = 1;  goto CLEANUP;
    }

    rhs = (double) c->rhs;
    coef = COLOR_SAFE_MALLOC (c->count, double);
    COLORcheck_NULL (coef, "out of memory for coef");
    for (i = 0; i < c->count; i++) coef[i] = 1.0;

    rval = COLORlp_addrow (lp, c->count, c->inodes, coef, COLORlp_LESS_EQUAL,
                           rhs, NULL);
    if (rval) COLORlp_printerrorcode (rval);
    COLORcheck_rval (rval, "COLORlp_addrow failed");

CLEANUP:

    COLOR_IFFREE (coef, double);
    return rval;
}

static void init_cut_pool (pool *P)
{
    P->cuts = (stablecut *) NULL;
    P->cutcount = 0;
    P->cutspace = 0;
}

static int add_cut_to_pool (pool *P, stablecut **pcset, int *hit)
{
    int rval = 0;

    if (P->cutcount >= P->cutspace) {
        printf ("Ran out of space in cut pool\n");
        rval = 1;  goto CLEANUP;
    }

    P->cuts[P->cutcount].inodes = (*pcset)->inodes;
    P->cuts[P->cutcount].count = (*pcset)->count;
    P->cuts[P->cutcount].type = (*pcset)->type;
    P->cuts[P->cutcount].rhs = (*pcset)->rhs;
    if (hit) *hit = P->cutcount;
    P->cutcount++;
    (*pcset)->inodes = (int *) NULL;

    COLOR_FREE (*pcset, stablecut);
    *pcset = (stablecut *) NULL;

CLEANUP:

    return rval;
}

static int find_violated_holes (COLORadjgraph *G, double *x, 
        stablecut **phlist, int *pcount)
{
    int i, j, k, ncount = G->ncount, rval = 0;
    int winner, trys, marker = 1, setcount = 0;
    int start, current, next, len, bestlen = 0, deg;
    int *nmarks = (int *) NULL;
    int *path = (int *) NULL;
    int *bestpath = (int *) NULL;
    int *ncover = (int *) NULL;
    COLORadjnode *n;
    double *nsum = (double *) NULL;
    double **eweights = (double **) NULL;
    double uval, t, pweight, bestviol, vmin = 1.0 - STABLE_EPS;
    stablecut *hcut = (stablecut *) NULL, *hlist = (stablecut *) NULL;
    COLORrandstate rstate;

    *phlist = (stablecut *) NULL;
    if (pcount) *pcount = 0;

    COLORutil_sprand (99, &rstate);

    eweights = COLOR_SAFE_MALLOC (ncount, double *);
    COLORcheck_NULL (eweights, "out of memory for eweights");
    for (i = 0; i < ncount; i++) eweights[i] = (double *) NULL;

    for (i = 0; i < ncount; i++) {
        eweights[i] = COLOR_SAFE_MALLOC (G->nodelist[i].degree, double);
        COLORcheck_NULL (eweights[i], "out of memory for eweights[i]");
    }

    nsum = COLOR_SAFE_MALLOC (ncount, double);
    COLORcheck_NULL (nsum, "out of memory for nsum");

    for (i = 0; i < ncount; i++) {
        nsum[i] = 0.0;
        n = &G->nodelist[i];
        for (j = 0; j < n->degree; j++) {
            eweights[i][j] = 1.0 - x[i] - x[n->adj[j]];
            nsum[i] += eweights[i][j];
        }
    }

    nmarks = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (nmarks, "out of memory for nmarks");
    for (i = 0; i < ncount; i++) nmarks[i] = 0;

    ncover = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (ncover, "out of memory for ncover");
    for (i = 0; i < ncount; i++) ncover[i] = 0;

    path = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (path, "out of memory for path");
    bestpath = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (bestpath, "out of memory for bestpath");


    for (start = 0; start < ncount; start++) {
        if (ncover[start]) continue;
        bestviol = 0.0;
        for (k = 0; k < 100; k++) {
            marker++;
            current = start;
            len = 0;
            pweight = 0.0;
            do {
                nmarks[current] = marker;
                path[len++] = current;
                winner = -1; trys = 0;
                deg = G->nodelist[current].degree;
                do {
                    trys++;
                    uval = ((double) COLORutil_lprand (&rstate)) / 
                           ((double) COLOR_PRANDMAX);
                    uval *= nsum[current];
                    t = 0.0;
                    for (j = 0; j < deg; j++) {
                        t += eweights[current][j];
                        if (t >= uval) {
                            next = G->nodelist[current].adj[j];
                            if ((nmarks[next] < marker ||
                                (next == start && len % 2 == 1 && len > 3)) &&
                                 pweight + eweights[current][j] < vmin) {
                                winner = next;
                                pweight += eweights[current][j];
                                break;
                            }
                        }
                    }
                } while (winner == -1 && trys < 2*deg);
                if (winner != start) current = winner;
            } while (winner != start && winner != -1);
            if (winner == start) {
                t = (1.0 - pweight) / (double) len;
                if (t > bestviol) {
                    bestviol = t;
                    for (i = 0; i < len; i++) {
                        bestpath[i] = path[i];
                    }
                    bestlen = len;
                }
            }
        }
        if (bestviol > 0.0) {
/*
            printf ("Best Odd-Hole for Node %d: Score = %.6lf\n",
                     start, bestviol);
            for (i = 0; i < bestlen; i++) printf ("%d ", bestpath[i]+1);
            printf ("\n"); fflush (stdout);
*/

            for (i = 0; i < bestlen; i++) ncover[bestpath[i]] = 1;
            rval = build_stablecut (bestpath, bestlen, STABLE_ODDHOLE, &hcut);
            COLORcheck_rval (rval, "build_stableset failed");
            hcut->next = hlist;
            hlist = hcut;
            setcount++;
            hcut = (stablecut *) NULL;
        }
    }

    *phlist = hlist;
    if (pcount) *pcount = setcount;


CLEANUP:

    COLOR_IFFREE (nsum, double);
    COLOR_IFFREE (ncover, int);
    COLOR_IFFREE (nmarks, int);
    COLOR_IFFREE (path, int);
    COLOR_IFFREE (bestpath, int);
    if (eweights) {
        for (i = 0; i < ncount; i++) {
            COLOR_IFFREE (eweights[i], double);
        }
        COLOR_FREE (eweights, double *);
    }
    return rval;
}

static void init_stable_prob (prob *p)
{
    p->icuts = (int *) NULL;
    p->cutcount = 0;
    p->cutspace = 0;
    p->branchdepth = 0;
    p->branchhistory = (branchobj *) NULL;
    p->warmstart = (COLORlp_warmstart *) NULL;
    p->next = (prob *) NULL;
    p->prev = (prob *) NULL;
} 

static void free_stable_prob (prob *p)
{
    if (p) {
        COLOR_IFFREE (p->icuts, int);
        COLOR_IFFREE (p->branchhistory, branchobj);
        if (p->warmstart) COLORlp_free_warmstart (&p->warmstart);
        init_stable_prob (p);
    }
}

static int build_stable_prob (prob *p, int setcount, int ncount)
{
    int rval = 0;

    p->cutcount = 0;
    p->cutspace = 2*setcount + 1000*ncount;
    p->icuts = COLOR_SAFE_MALLOC (p->cutspace, int);
    COLORcheck_NULL (p->icuts, "out of memory for icuts");

CLEANUP:

    return rval;
}

static int add_cut_to_prob (prob *p, int ind)
{
    int rval = 0;

    if (p->cutcount >= p->cutspace) {
        printf ("Ran out of space in cutlist\n");
        rval = 1;  goto CLEANUP;
    }
    p->icuts[p->cutcount] = ind;
    p->cutcount++;

CLEANUP:

    return rval;
}

static int build_cut_pool (pool *P, int setcount, int ncount)
{
    int i, rval = 0;

    P->cutspace = 2*setcount + 1000*ncount;
    P->cuts = COLOR_SAFE_MALLOC (P->cutspace, stablecut);
    COLORcheck_NULL (P->cuts, "out of memory for cuts");
    for (i = 0; i < P->cutspace; i++) init_stablecut (&P->cuts[i]);
    P->cutcount = 0;

CLEANUP:

    return rval;
}

static void free_cut_pool (pool *P)
{
    int i;

    for (i = 0; i < P->cutcount; i++) {
        free_stablecut (&P->cuts[i]);
    }
    COLOR_IFFREE (P->cuts, stablecut);
}

static int cover_edge_cliques (COLORadjgraph *G, int ecount, int *elist,
        stablecut **pclist, int *setcount)
{
    int rval = 0;
    int i, count = 0, ncount = G->ncount;
    int **incid = (int **) NULL;
    int *emarks = (int *) NULL;
    int *nmarks = (int *) NULL;
    int *nlist = (int *) NULL;
    stablecut *cset = (stablecut *) NULL;
    stablecut *clist = (stablecut *) NULL;

    *pclist = (stablecut *) NULL;
    *setcount = 0;

    incid = COLOR_SAFE_MALLOC (ncount, int *);
    COLORcheck_NULL (incid, "out of memory for incid"); 

    for (i = 0; i < ncount; i++) incid[i] = (int *) NULL;
    rval = build_incidence (G, ecount, elist, incid);
    COLORcheck_rval (rval, "build_incidence failed");

    emarks = COLOR_SAFE_MALLOC (ecount, int);
    COLORcheck_NULL (emarks, "out of memory for emarks");

    nmarks = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (nmarks, "out of memory for nmarks"); 

    nlist = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (nlist, "out of memory for nlist"); 

    for (i = 0; i < ecount; i++) emarks[i] = 0;
    for (i = 0; i < ncount; i++) nmarks[i] = 0;

    for (i = 0; i < ecount; i++) {
        if (emarks[i] == 0) {
            rval = cover_edge (G, elist[2*i], elist[2*i+1], &cset, nmarks,
                               nlist);
            COLORcheck_rval (rval, "cover_edge failed");
            cset->next = clist;
            clist = cset;
            count++;
            mark_edges (G, incid, cset, emarks, nmarks, elist); 
        }
    } 

    *pclist = clist;
    *setcount = count;

CLEANUP:

    COLOR_IFFREE (emarks, int);
    COLOR_IFFREE (nmarks, int);
    COLOR_IFFREE (nlist, int);
    if (incid) {
        for (i = 0; i < ncount; i++) {
            COLOR_IFFREE (incid[i], int);
        }
        COLOR_FREE (incid, int *);
    }

    return rval;
}

static int find_violated_cliques (COLORadjgraph *G, double *x, 
        stablecut **pclist, int *pcount)
{
    int rval = 0;
    int i, j, k, yesno, count = 0, ncount = G->ncount;
    int *nlist = (int *) NULL;
    int *nmarks = (int *) NULL;
    int *cover = (int *) NULL;
    int *perm = (int *) NULL;
    stablecut *cset = (stablecut *) NULL;
    stablecut *clist = (stablecut *) NULL;
    double val;

    *pclist = (stablecut *) NULL;

    nlist = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (nlist, "out of memory for nlist");

    nmarks = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (nmarks, "out of memory for nmarks");
    for (i = 0; i < ncount; i++) nmarks[i] = 0;

    cover = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (cover, "out of memory for cover");
    for (i = 0; i < ncount; i++) cover[i] = 0;

    perm = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (perm, "out of memory for perm");
    for (i = 0; i < ncount; i++) perm[i] = i;
    perm_dbl_rquicksort (perm, x, ncount);

    for (i = 0; i < ncount && x[perm[i]] > STABLE_EPS; i++) {
        k = perm[i];
        if (cover[k] == 0) {
            rval = cover_node (G, x, k, &cset, nmarks, nlist);
            COLORcheck_rval (rval, "cover_node failed");
            if (cset) {
                val = 0;
                for (j = 0; j < cset->count; j++) {
                    val += x[cset->inodes[j]];
                }
                if (val > 1.0 + STABLE_VIOL) {
                    /* printf ("Clique %d: %.6lf\n", cset->count, val); */
                    check_clique (G, nmarks, cset, &yesno);
                    if (yesno == 1) {
                        cset->next = clist;
                        clist = cset;
                        for (j = 0; j < cset->count; j++) {
                            cover[cset->inodes[j]] = 1;
                        }
                        cset = (stablecut *) NULL;
                        count++;
                    }
                }
                if (cset) {
                    free_stablecut (cset);
                    COLOR_FREE (cset, stablecut);
                }
            }
        }
    }

    *pclist = clist;
    if (pcount) *pcount = count;

CLEANUP:

    COLOR_IFFREE (nlist, int);
    COLOR_IFFREE (nmarks, int);
    COLOR_IFFREE (cover, int);
    COLOR_IFFREE (perm, int);
    return rval;
}

static int cover_node (COLORadjgraph *G, double *weights, int k,
        stablecut **pcset, int *nmarks, int *nlist)
{
    int rval = 0;
    COLORadjnode *n;
    int i, winner, count;
    double bestweight;

    *pcset = (stablecut *) NULL;

    count = 0;
    winner = k;
    while (winner != -1) {
        n = &G->nodelist[winner];
        nlist[count] = winner;
        count++;
        winner = -1;
        bestweight = 0.0;
        for (i = 0; i < n->degree; i++) {
            if (nmarks[n->adj[i]] == count-1) {
                if (weights[n->adj[i]] > bestweight) {
                    bestweight = weights[n->adj[i]];
                    winner = n->adj[i];
                }
                nmarks[n->adj[i]] = count;
            }
        }
    }

    if (count >= 3) {
        rval = build_stablecut (nlist, count, STABLE_CLIQUE, pcset);
        COLORcheck_rval (rval, "build_stablecut failed");
    }

CLEANUP:

    for (i = 0; i < G->ncount; i++) nmarks[i] = 0;
    return rval;
}

static int cover_edge (COLORadjgraph *G, int end1, int end2, stablecut **pcset,
        int *nmarks, int *nlist)
{
    int rval = 0;
    COLORadjnode *n;
    int i, count, winner;

    nlist[0] = end1;
    n = &G->nodelist[end1];
    for (i = 0; i < n->degree; i++) {
        nmarks[n->adj[i]] = 1;
    }

    count = 1;
    winner = end2;
    while (winner != -1) {
        n = &G->nodelist[winner];
        nlist[count] = winner;
        count++;
        winner = -1;
        for (i = 0; i < n->degree; i++) {
            if (nmarks[n->adj[i]] == count-1) {
                winner = n->adj[i];
                nmarks[winner] = count;
            }
        }
    }

    rval = build_stablecut (nlist, count, STABLE_CLIQUE, pcset);
    COLORcheck_rval (rval, "build_stable failed");

CLEANUP:

    for (i = 0; i < G->ncount; i++) nmarks[i] = 0;
    return rval;
}

static void check_clique (COLORadjgraph *G, int *nmarks, stablecut *cset, 
        int *yesno)
{
    int i, j;
    COLORadjnode *n;

    *yesno = 0;
    for (i = 0; i < G->ncount; i++) nmarks[i] = 0;

    for (i = 0; i < cset->count; i++) {
        if (nmarks[cset->inodes[i]]) {
            printf ("Duplicate node in clique\n");
            goto CLEANUP;
        }
        nmarks[cset->inodes[i]] = 1;
    }
    for (i = 0; i < cset->count; i++) nmarks[cset->inodes[i]] = 0;

    for (i = 0; i < cset->count; i++) {
        n = &G->nodelist[cset->inodes[i]];
        for (j = 0; j < n->degree; j++) {
            nmarks[n->adj[j]]++;
        }
    }

    for (i = 0; i < cset->count; i++) {
        if (nmarks[cset->inodes[i]] != cset->count - 1) {
            printf ("Missing degree in clique\n");
            goto CLEANUP;
        }
    }

    *yesno = 1;

CLEANUP:

    for (i = 0; i < G->ncount; i++) nmarks[i] = 0;
    if (*yesno == 0) {
        printf ("Set is not a clique\n");
    }
}

static void mark_edges (COLORadjgraph *G, int **incid, stablecut *cset, 
        int *emarks, int *nmarks, int *elist) 
{
    int i, j, k, n, n1, n2;

    for (i = 0; i < cset->count; i++) nmarks[cset->inodes[i]] = 1;
    for (i = 0; i < cset->count; i++) {
        n = cset->inodes[i];
        for (j = 0; j < G->nodelist[n].degree; j++) {
            k = incid[n][j];
            n1 = elist[2*k];
            n2 = elist[2*k+1];
            if (nmarks[n1] == 1 && nmarks[n2] == 1) {
                emarks[k] = 1;
            }
        }
    }
    for (i = 0; i < cset->count; i++) nmarks[cset->inodes[i]] = 0;
}

static int stable_greedy_x (COLORadjgraph *G, int *weights, double *x,
        COLORset **pcset, int *wt)
{
    int rval = 0;
    int i, j, w, w2, count = 0, count2 = 0, ncount = G->ncount;
    int *perm = (int *) NULL;
    int *nmarks = (int *) NULL;
    int *nlist = (int *) NULL, *nlist2 = (int *) NULL;
    double *nx = (double *) NULL;
    COLORadjnode *n;

    perm = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (perm, "out of memory for perm");
    for (i = 0; i < ncount; i++) perm[i] = i;
    perm_dbl_rquicksort (perm, x, ncount);

    nmarks = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (nmarks, "out of memory for nmarks");
    for (i = 0; i < ncount; i++) nmarks[i] = 0;

    nlist = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (nlist, "out of memory for nlist");
    nlist2 = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (nlist, "out of memory for nlist2");

    for (i = 0; i < ncount; i++) {
        if (nmarks[perm[i]] == 0) {
            nmarks[perm[i]] = 1;
            nlist[count++] = perm[i];
            n = &G->nodelist[perm[i]];
            for (j = 0; j < n->degree; j++) {
                nmarks[n->adj[j]] = 1;
            }
        }
    }

    w = 0;
    for (i = 0; i < count; i++) w += weights[nlist[i]];

    nx = COLOR_SAFE_MALLOC (ncount, double);
    COLORcheck_NULL (nx, "out of memory for nx");
    for (i = 0; i < ncount; i++) nx[i] = x[i];
    for (i = 0; i < ncount; i++) {
        n = &G->nodelist[i];
        for (j = 0; j < n->degree; j++) {
            nx[i] += x[n->adj[j]];
        }
    }
    for (i = 0; i < ncount; i++) {
        nx[i] *= -1.0;
        perm[i] = i;
    }
    perm_dbl_rquicksort (perm, nx, ncount);
    for (i = 0; i < ncount; i++) nmarks[i] = 0;

    for (i = 0; i < ncount; i++) {
        if (nmarks[perm[i]] == 0) {
            nmarks[perm[i]] = 1;
            nlist2[count2++] = perm[i];
            n = &G->nodelist[perm[i]];
            for (j = 0; j < n->degree; j++) {
                nmarks[n->adj[j]] = 1;
            }
        }
    }

    w2 = 0;
    for (i = 0; i < count2; i++) w2 += weights[nlist2[i]];

    if (w > w2) {
        *wt = w;
        rval = build_set (nlist, count, pcset);
    } else {
        *wt = w2;
        rval = build_set (nlist2, count2, pcset);
    }
    COLORcheck_rval (rval, "build_set failed");

CLEANUP:

    COLOR_IFFREE (perm, int);
    COLOR_IFFREE (nmarks, int);
    COLOR_IFFREE (nlist, int);
    COLOR_IFFREE (nlist2, int);
    COLOR_IFFREE (nx, double);
    return rval;
}

static int call_dfs_branching (COLORadjgraph *G, int *weights, COLORlp *lp,
        pool *P, prob *probdata, int *bestval)
{
    int rval = 0;
    int bcount = 0;
    double *x = (double *) NULL;

    printf ("Call DFS Branching ...\n"); fflush (stdout);

    x = COLOR_SAFE_MALLOC (G->ncount, double);
    COLORcheck_NULL (x, "out of memory for x");

    rval = dfs_branching (G, weights, lp, P, probdata, bestval, x, 0, &bcount);
    COLORcheck_rval (rval, "dfs_branching failed");

    printf ("Best Value: %d\n", *bestval);

CLEANUP:

    COLOR_IFFREE (x, double);
    return rval;
}

static int dfs_branching (COLORadjgraph *G, int *weights, COLORlp *lp, pool *P,
        prob *probdata, int *bestval, double *x, int depth, int *bcount)
{
    int rval = 0;
    int bvertex;
    double val;

    (*bcount)++;

    rval = call_lp_solver (lp, &val, x);
    COLORcheck_rval (rval, "call_lp_solver failed");

    printf ("Depth %2d, BBnodes %d, LP Value: %.6f ", depth, *bcount, val);

    if (val - ((double) *bestval) < (1.0 - STABLE_EPS)) {
        printf ("Prune\n");
        goto CLEANUP;
    }

    rval = find_branch (G, x, &bvertex);
    COLORcheck_rval (rval, "find_branch failed");

    if (bvertex == -1) {
        COLORset cset;
        int solval;

        COLORinit_set (&cset);
        rval = get_integer_soln (G, weights, x, &cset, &solval);
        COLORcheck_rval (rval, "get_integer_soln failed");

        if (solval > *bestval) {
            *bestval = solval;
            printf ("Integer Soln: %d\n", solval); fflush (stdout);
        }

        COLORfree_set (&cset);
        goto CLEANUP;
    }

    printf ("\n");  fflush (stdout);

    /* 1-Side Branch */

    rval = COLORlp_setbound (lp, bvertex, 'L', 1.0);
    COLORcheck_rval (rval, "set_branch failed");

    rval = cutting_loop (G, lp, P, probdata, 1);
    COLORcheck_rval (rval, "cutting_loop failed");

    rval = dfs_branching (G, weights, lp, P, probdata, bestval, x, depth+1,
                          bcount);
    COLORcheck_rval (rval, "dfs_branching failed");

    /* 0-Side Branch */

    rval = COLORlp_setbound (lp, bvertex, 'L', 0.0);
    COLORcheck_rval (rval, "set_branch failed");

    rval = COLORlp_setbound (lp, bvertex, 'U', 0.0);
    COLORcheck_rval (rval, "set_branch failed");

    rval = cutting_loop (G, lp, P, probdata, 1);
    COLORcheck_rval (rval, "cutting_loop failed");

    rval = dfs_branching (G, weights, lp, P, probdata, bestval, x, depth+1,
                          bcount);
    COLORcheck_rval (rval, "dfs_branching failed");

    rval = COLORlp_setbound (lp, bvertex, 'U', 1.0);
    COLORcheck_rval (rval, "set_branch failed");

CLEANUP:

    return rval;
}

static int find_branch (COLORadjgraph *G, double *x, int *pvertex)
{
    int rval = 0;
    int i, ncount = G->ncount;
    double *dist = (double *) NULL;
    int *perm = (int *) NULL;

    *pvertex = -2;

    dist = COLOR_SAFE_MALLOC (ncount, double);
    COLORcheck_NULL (dist, "out of memory for dist");

    for (i = 0; i < ncount; i++) {
        if (x[i] > 0.5) dist[i] = x[i] - 0.5;
        else            dist[i] = 0.5 - x[i];
    }

    perm = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (perm, "out of memory for perm");
    for (i = 0; i < ncount; i++) perm[i] = i;
    perm_dbl_quicksort (perm, dist, ncount);

    if (dist[perm[0]] > 0.5 - STABLE_EPS) {
        /*  LP solution is integer valued */
        *pvertex = -1;
    } else {
        *pvertex = perm[0];
    }
    
CLEANUP:

    COLOR_IFFREE (dist, double);
    COLOR_IFFREE (perm, int);
    return rval;
}

static int get_integer_soln (COLORadjgraph *G, int *weights, double *x, 
        COLORset *sol, int *solval)
{
    int i, w = 0, cnt = 0, rval = 0;

    *solval = 0;

    for (i = 0; i < G->ncount; i++) {
        if (x[i] > 1.0 - STABLE_EPS) cnt++;
    }

    sol->members = COLOR_SAFE_MALLOC (cnt, int);
    COLORcheck_NULL (sol->members, "out of memory for members");

    sol->count = 0;
    for (i = 0; i < G->ncount; i++) {
        if (x[i] > 1.0 - STABLE_EPS) {
            sol->members[sol->count++] = i;
            w += weights[i];
        }
    }

    *solval = w;

CLEANUP:

    return rval;
}

static int build_incidence (COLORadjgraph *G, int ecount, int *elist, 
        int **incid)
{
    int rval = 0;
    int i, n1, n2;
    int *deg = (int *) NULL;

    for (i = 0; i < G->ncount; i++) {
        incid[i] = COLOR_SAFE_MALLOC (G->nodelist[i].degree, int);
        COLORcheck_NULL (incid[i], "out of memory for incid");
    } 

    deg = COLOR_SAFE_MALLOC (G->ncount, int);
    COLORcheck_NULL (deg, "out of memory for deg");
    for (i = 0; i < G->ncount; i++) deg[i] = 0;

    for (i = 0; i < ecount; i++) {
        n1 = elist[2*i];
        n2 = elist[2*i+1];
        incid[n1][deg[n1]++] = i;
        incid[n2][deg[n2]++] = i;
    }

CLEANUP:

    COLOR_IFFREE (deg, int);
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

   return rval;
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

