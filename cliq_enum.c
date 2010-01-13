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
**/


/****************************************************************************/
/*                                                                          */
/*  EXPORTED FUNCTIONS:                                                     */
/*                                                                          */
/*  int COLORclique_ostergard (COLORset **newsets, int *nnewsets,           */
/*      int ncount, int ecount, int *elist, int *weights, int cutoff,       */
/*      int *pval)                                                          */
/*    FINDS a maximum-weight clique via P. Ostergard's enumeration alg.     */
/*      The code is really only suitable for sparse graphs.                 */
/*    -newsets returns the clique (or cliques, but only one currently)      */
/*    -nnewsets returns the number of cliques (this will be one currently)  */
/*    -ncount is the number of nodes in the graph                           */
/*    -ecount is the number of edges in the graph                           */
/*    -elist is an array of the edges in end0 end1 format (2*ecount long)   */
/*    -weights is an array of node weighs (ncount long)                     */
/*    -cutoff specifies a target weight for the clique; the code will stop  */
/*     the search once a clique of weight cutoff is found; to find an       */
/*     optimal clique set cutoff to COLOR_MAXINT                            */
/*    -pval returns the weight of the clique (can be NULL)                  */
/*                                                                          */
/*  int COLORclique_enum (COLORset** newsets, int *nnewsets, int ncount,    */
/*      int ecount, int *elist, int *weights, int cutoff, int *pval);       */
/*    FINDS a maximum-weight weight clique by enumeration.  The method is   */
/*      only suitable for sparse graphs.                                    */
/*    -parameters same as COLORclique_ostergard                             */
/*                                                                          */
/****************************************************************************/


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

static int run_clique_enum (int ncount, int ecount, int *elist, int *weights,
    int *pval, int first, int *marks, int *invmap, int cutoff,
    COLORset *bestcliq);
static int main_ostergard (int ncount, int ecount, int *elist, int *weights,
        int *optval, int *bestcnt, int *bestset, int cutoff);
static int run_ostergard (COLORadjgraph *G, int Ucount, int *U, int *weights, int psum,
        int *C, int *bigmax, int *marks, int *pcnt, int *pset, int *bestcnt,
        int *bestset);
static int grab_marked_neighbors (COLORadjgraph *G, int *marks, int k, int *pcount,
        int **plist);
static int ostergard_order (COLORadjgraph *G, int *weights, int *order);
static int permute_nodes (int *invorder, int ncount, int ecount, int *elist,
        int *weights, int **pielist, int **piweights);
static int build_directed_adj (COLORadjgraph * G, int ncount, int ecount, int *elist,
        int forward);
static int check_clique (int ncount, int ecount, int *elist, COLORset *clique,
    int *yesno);
static void max_unmarked (int ncount, int *marks, int *weights, int *sweights, 
        int *k);

int COLORclique_enum (COLORset **newsets, int *nnewsets, int ncount,
        int ecount, int *elist, int *weights, int cutoff, int *pval)
{
    int rval = 0;
    int i, yesno = 0, val = 0;
    int *marks = (int *) NULL, *invmap = (int *) NULL;
    COLORset cliq;

    if (nnewsets) *nnewsets = 0;
    if (newsets) *newsets = (COLORset *) NULL;

    COLORinit_set (&cliq);
    cliq.members = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (cliq.members, "out of memory for cliq.members");

    for (i = 0; i < ecount; i++) {
        if (elist[2*i] > elist[2*i+1]) {
            printf ("Error: reversed edge in elist\n"); fflush (stdout);
            rval = 1; goto CLEANUP;
        }
        if (i > 1) {
            if (elist[2*i] < elist[2*i-2] || (elist[2*i] == elist[2*i-2] && 
                elist[2*i+1] < elist[2*i-1])) {
                printf ("Error: elist not sorted\n"); fflush (stdout);
                rval = 1; goto CLEANUP;
            }
        }
    }

    marks = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (marks, "out of memory for marks");
    for (i = 0; i < ncount; i++) marks[i] = 0;

    invmap = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (marks, "out of memory for invmap");

    rval = run_clique_enum (ncount, ecount, elist, weights, &val, 0, marks,
                            invmap, cutoff, &cliq);
    COLORcheck_rval (rval, "run_clique_enum failed");

    rval = check_clique (ncount, ecount, elist, &cliq, &yesno);
    COLORcheck_rval (rval, "check_clique failed");
    
    if (yesno == 0) {
        printf ("ERROR: run_clique_enum returned without a clique\n");
        rval = 1;  goto CLEANUP;
    }

    if (pval) *pval = val;

    *newsets = COLOR_SAFE_MALLOC (1, COLORset);
    COLORcheck_NULL (*newsets, "out of memory for newsets");
    COLORinit_set (*newsets);
    (*newsets)->count = cliq.count;
    (*newsets)->members = cliq.members;
    *nnewsets = 1;
    cliq.members = (int *) NULL;

CLEANUP:

    COLOR_IFFREE (marks, int);
    COLOR_IFFREE (invmap, int);
    COLORfree_set (&cliq);
    return rval;
}

static int run_clique_enum (int ncount, int ecount, int *elist, int *weights,
        int *pval, int first, int *marks, int *invmap, int cutoff,
        COLORset *bestcliq)
{
    int rval = 0;
    int i, hval = 0, nval = 0, deg, hecount;
    int bestval = -1, best = -1;
    int *pe = elist;
    int *hweights = (int *) NULL, *helist = (int *) NULL;
    int *map = (int *) NULL;
    COLORset hcliq;

    COLORinit_set (&hcliq);

    if (ncount == 1) {
        *pval = weights[first];
        bestcliq->members[0] = first;
        bestcliq->count = 1;      
        goto CLEANUP;
    }

    if (ecount == 0) {
        for (i = 0; i < ncount; i++) {
            if (weights[i+first] > bestval) {
                bestval = weights[i+first];
                best = i+first;
            }
        }
        *pval = bestval;
        bestcliq->members[0] = best;
        bestcliq->count = 1;      
        goto CLEANUP;
    }

    hcliq.members = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (hcliq.members, "out of memory for hcliq.members");

    nval = weights[first];
    deg = 0;
    for (i = 0; i < ecount; i++) {
        if (elist[2*i] != first) break;
        deg++;
    }

    if (deg) {
        hweights = COLOR_SAFE_MALLOC (deg, int);
        COLORcheck_NULL (hweights, "out of memory for hweights");

        map = COLOR_SAFE_MALLOC (deg, int);
        COLORcheck_NULL (hweights, "out of memory for map");

        for (i = 0; i < deg; i++) {
            map[i] = elist[2*i+1];
            hweights[i] = weights[map[i]];
            marks[map[i]] = 1;
            invmap[map[i]] = i;
        }

        hecount = 0;
        for (i = deg; i < ecount; i++) {
            if (marks[elist[2*i]] == 1 && marks[elist[2*i+1]] == 1) {
                hecount++;
            }
        }

        if (hecount) {
            helist = COLOR_SAFE_MALLOC (2*hecount, int);
            COLORcheck_NULL (helist, "out of memory for helist");

            hecount = 0;
            for (i = deg; i < ecount; i++) {
                if (marks[elist[2*i]] == 1 && marks[elist[2*i+1]] == 1) {
                    helist[2*hecount] = invmap[elist[2*i]];
                    helist[2*hecount+1] = invmap[elist[2*i+1]];
                    hecount++;
                }
            }
        }

        for (i = 0; i < deg; i++) marks[map[i]] = 0;

        rval = run_clique_enum (deg, hecount, helist, hweights, &hval, 0,
                                marks, invmap, cutoff, &hcliq);
        COLORcheck_rval (rval, "run_clique_enum failed");
 
        if (nval + hval > bestval) {
            bestval = nval + hval;
            bestcliq->count = hcliq.count + 1;      
            for (i = 0; i < hcliq.count; i++) {
                bestcliq->members[i] = map[hcliq.members[i]];
            }
            bestcliq->members[i] = first;
        }
    } else {
        if (nval > bestval) {
            bestval = nval;
            bestcliq->members[0] = first;
            bestcliq->count = 1;      
        }
    }

    if (bestval > cutoff) {
        *pval = bestval;
        goto CLEANUP;
    }

    ncount--;
    ecount -= deg;
    pe += (2*deg);
      
    rval = run_clique_enum (ncount, ecount, pe, weights, &hval, first+1,
                            marks, invmap, cutoff, &hcliq);
    COLORcheck_rval (rval, "run_clique_enum failed");

    if (hval > bestval) {
        bestval = hval;
        bestcliq->count = hcliq.count;      
        for (i = 0; i < hcliq.count; i++) {
            bestcliq->members[i] = hcliq.members[i];
        }
    }

    *pval = bestval;

CLEANUP:

    COLOR_IFFREE (hweights, int);
    COLOR_IFFREE (helist, int);
    COLOR_IFFREE (map, int);
    COLORfree_set (&hcliq);
    return rval;
}

int COLORclique_ostergard (COLORset **newsets, int *nnewsets, int ncount,
        int ecount, int *elist, int *weights, int cutoff, int *pval)
{
    int i, yesno = 0, rval = 0, optval = 0;
    int *order = (int *) NULL, *invorder = (int *) NULL;
    int *ielist = (int *) NULL, *iweights = (int *) NULL;
    int bestcnt = 0, *bestset = (int *) NULL;
    COLORadjgraph G;

    COLORadjgraph_init (&G);

    if (nnewsets) *nnewsets = 0;
    if (newsets) *newsets = (COLORset *) NULL;
    if (pval) *pval = 0;

    rval = COLORadjgraph_build (&G, ncount, ecount, elist);
    COLORcheck_rval (rval, "COLORadjgraph_build failed");

    order = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (order, "out of memory for order");

    rval = ostergard_order (&G, weights, order);
    COLORcheck_rval (rval, "ostergard_order failed");

    invorder = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (invorder, "out of memory for invorder");
    for (i = 0; i < ncount; i++) invorder[order[i]] = i;

    bestset = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (bestset, "out of memory for bestset");

    rval = permute_nodes (invorder, ncount, ecount, elist, weights,
                          &ielist, &iweights);
    COLORcheck_rval (rval, "permute nodes failed");

    rval = main_ostergard (ncount, ecount, ielist, iweights, &optval,
                           &bestcnt, bestset, cutoff);
    COLORcheck_rval (rval, "main_ostergard failed");

    *newsets = COLOR_SAFE_MALLOC (1, COLORset);
    COLORcheck_NULL (*newsets, "out of memory for newsets");
    COLORinit_set (*newsets);
    (*newsets)->count = bestcnt;
    (*newsets)->members = COLOR_SAFE_MALLOC (bestcnt, int);
    for (i = 0; i < bestcnt; i++) {
        (*newsets)->members[i] = order[bestset[i]];
    }
    COLORutil_quicksort ((*newsets)->members, bestcnt);
    *nnewsets = 1;

    if (pval) *pval = optval;

    rval = check_clique (ncount, ecount, elist, *newsets, &yesno);
    COLORcheck_rval (rval, "check_clique failed");
    
    if (yesno == 0) {
        printf ("ERROR: ostergard code returned without a clique\n");
        rval = 1;  goto CLEANUP;
    }

CLEANUP:

    COLOR_IFFREE (order, int);
    COLOR_IFFREE (invorder, int);
    COLOR_IFFREE (ielist, int);
    COLOR_IFFREE (iweights, int);
    COLOR_IFFREE (bestset, int);
    COLORadjgraph_free (&G);
    return rval;
}

static int main_ostergard (int ncount, int ecount, int *elist, int *weights,
        int *optval, int *bestcnt, int *bestset, int cutoff)
{
    int i, rval = 0;
    int *C = (int *) NULL, *marks = (int *) NULL;
    int bigmax = -COLOR_MAXINT, *U, Ucount;
    int pcnt = 0, *pset = (int *) NULL;
    COLORadjgraph G;

    COLORadjgraph_init (&G);

    pset = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (pset, "out of memory for pset");

    rval = build_directed_adj (&G, ncount, ecount, elist, 0);
    COLORcheck_rval (rval, "build_directed_adj failed");

    marks = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (marks, "out of memory for marks");
    for (i = 0; i < ncount; i++) marks[i] = 0;

    C = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (C, "out of memory for C");

    C[0] = weights[0]; bigmax = weights[0];
    *bestcnt = 1; bestset[0] = 0;
    if (optval) *optval = bigmax;
    if (bigmax >= cutoff) goto CLEANUP;

    for (i = 1; i < ncount && bigmax < cutoff; i++) {
        if (G.nodelist[i].degree == 0) {
            if (weights[i] > bigmax) {
                bigmax = weights[i];
            }
            C[i] = bigmax;
        } else {
            Ucount = G.nodelist[i].degree;
            U = G.nodelist[i].adj;

            pset[pcnt++] = i;
            run_ostergard (&G, Ucount, U, weights, weights[i], C, &bigmax,
                           marks, &pcnt, pset, bestcnt, bestset);
            COLORcheck_rval (rval, "run_ostergard failed");
            pcnt--;
            C[i] = bigmax;
        }
    }

    if (optval) *optval = bigmax;

CLEANUP:  

    COLORadjgraph_free (&G);
    COLOR_IFFREE (C, int);
    COLOR_IFFREE (marks, int);
    COLOR_IFFREE (pset, int);
    return rval;
}

static int run_ostergard (COLORadjgraph *G, int Ucount, int *U, int *weights, int psum,
        int *C, int *bigmax, int *marks, int *pcnt, int *pset, int *bestcnt,
        int *bestset)
{
    int rval = 0;
    int i, j, k, m, tsum = 0;
    int newUcount = 0, *newU = (int *) NULL;

    if (Ucount == 1) {
        if (psum + weights[U[0]] > *bigmax) {
            *bigmax = psum+weights[U[0]];
            for (m = 0; m < *pcnt; m++) {
                bestset[m] = pset[m];
            }
            bestset[m] = U[0];
            *bestcnt = *pcnt + 1;
        }
        goto CLEANUP;
    }

    for (i = 0, tsum = 0; i < Ucount; i++) tsum += weights[U[i]];

    for (i = 0; i < Ucount; i++) {
        if (psum + tsum <= *bigmax || psum + C[U[i]] <= *bigmax) break;

        k = U[i];

        for (j = i+1; j < Ucount; j++) marks[U[j]] = 1;
        rval = grab_marked_neighbors (G, marks, k, &newUcount, &newU);
        COLORcheck_rval (rval, "grab_marked_neighbors failed");
        for (j = i+1; j < Ucount; j++) marks[U[j]] = 0;

        if (newUcount) {
            pset[(*pcnt)++] = k;
            rval = run_ostergard (G, newUcount, newU, weights,
                         psum+weights[k], C, bigmax, marks, pcnt, pset,
                         bestcnt, bestset);
            COLORcheck_rval (rval, "run_ostergard failed"); 
            (*pcnt)--;
        } else {
            if (psum + weights[k] > *bigmax) {
                *bigmax = psum+weights[k];
                for (m = 0; m < *pcnt; m++) {
                    bestset[m] = pset[m];
                }
                bestset[m] = k;
                *bestcnt = *pcnt + 1;
            }
        }

        COLOR_IFFREE (newU, int);
        tsum -= weights[k];
    }

CLEANUP:

    COLOR_IFFREE (newU, int);
    return rval;
}

static int grab_marked_neighbors (COLORadjgraph *G, int *marks, int k, int *pcount,
        int **plist)
{
    int i, rval = 0;
    COLORadjnode *n = &G->nodelist[k];
    int count = 0, *list = (int *) NULL;

    *pcount = 0;
    *plist = (int *) NULL;

    if (n->degree == 0) goto CLEANUP;

    list = COLOR_SAFE_MALLOC (n->degree, int);
    COLORcheck_NULL (list, "out of memory for list");

    for (i = 0; i < n->degree; i++) {
        if (marks[n->adj[i]] == 1) {
            list[count++] = n->adj[i];
        }
    } 
    if (count == 0) {
        COLOR_FREE (list, int);
        goto CLEANUP;
    }

    *plist = list;
    *pcount = count;

CLEANUP:

    if (rval) {
        COLOR_IFFREE (list, int);
    }
    return rval;
}

static int ostergard_order (COLORadjgraph *G, int *weights, int *order)
{
    int rval = 0;
    int i, j, k, *sweights = (int *) NULL, ncount = G->ncount;
    int *marks = (int *) NULL;
    COLORadjnode *n;

    sweights = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (sweights, "out of memory for sweights");

    for (i = 0; i < ncount; i++) {
        sweights[i] = 0;
        n = &G->nodelist[i];
        for (j = 0; j < n->degree; j++) {
            sweights[i] += weights[n->adj[j]];
        }
    }

    marks = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (marks, "out of memory for marks");
    for (i = 0; i < ncount; i++) marks[i] = 0;

    for (i = 0; i < ncount; i++) {
        max_unmarked (ncount, marks, weights, sweights, &k);
        marks[k] = 1;
        order[i] = k;
        n = &G->nodelist[k];
        for (j = 0; j < n->degree; j++) {
            sweights[n->adj[j]] -= weights[i];
        }
    }

CLEANUP:

    COLOR_IFFREE (sweights, int);
    COLOR_IFFREE (marks, int);
    return rval;
}

static int permute_nodes (int *invorder, int ncount, int ecount, int *elist,
        int *weights, int **pielist, int **piweights)
{
    int i, rval = 0;
    int *ielist = (int *) NULL, *iweights = (int *) NULL;

    *pielist = (int *) NULL;
    *piweights = (int *) NULL;

    ielist = COLOR_SAFE_MALLOC (2*ecount, int);
    COLORcheck_NULL (pielist, "out of memory for pielist");

    iweights = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (iweights, "out of memory for iweights");

    for (i = 0; i < ecount; i++) {
        if (invorder[elist[2*i]] < invorder[elist[2*i+1]]) {
            ielist[2*i] = invorder[elist[2*i]];
            ielist[2*i+1] = invorder[elist[2*i+1]];
        } else {
            ielist[2*i] = invorder[elist[2*i+1]];
            ielist[2*i+1] = invorder[elist[2*i]];
        }
    }

    for (i = 0; i < ncount; i++) {
        iweights[invorder[i]] = weights[i];
    }

    *pielist = ielist;
    *piweights = iweights;

CLEANUP:

    if (rval) {
        COLOR_IFFREE (ielist, int);
        COLOR_IFFREE (iweights, int);
    }
    return rval;
}

static void max_unmarked (int ncount, int *marks, int *weights, int *sweights, 
        int *k)
{
    int i, wmin = COLOR_MAXINT, wmax = -COLOR_MAXINT;

    *k = -1;
    
    for (i = 0; i < ncount; i++) {
        if (marks[i] == 0) {
            if (weights[i] < wmin || (weights[i] == wmin &&
                                     sweights[i] > wmax)) {
                wmin = weights[i];
                wmax = sweights[i];
                *k = i;
            }
        }
    }
}

static int check_clique (int ncount, int ecount, int *elist, COLORset *clique,
        int *yesno)
{
    int i, rval = 0;
    int *nmarks = (int *) NULL;

    *yesno = 0;

    nmarks = COLOR_SAFE_MALLOC (ncount, int);
    COLORcheck_NULL (nmarks, "out of memory for nmarks");
    for (i = 0; i < ncount; i++) nmarks[i] = 0;

    for (i = 0; i < clique->count; i++) {
        if (nmarks[clique->members[i]]) {
            printf ("Duplicate node in clique\n");
            goto CLEANUP;
        }
        nmarks[clique->members[i]] = 1;
    }

    for (i = 0; i < ecount; i++) {
        if (nmarks[elist[2*i]] && nmarks[elist[2*i+1]]) {
            nmarks[elist[2*i]]++;
            nmarks[elist[2*i+1]]++;
        }
    }
     
    for (i = 0; i < clique->count; i++) {
        if (nmarks[clique->members[i]] != clique->count) {
            printf ("Missing degree in clique\n");
            goto CLEANUP;
        }
    }

    *yesno = 1;

CLEANUP:

    COLOR_IFFREE (nmarks, int);
    return rval;
}

static int build_directed_adj (COLORadjgraph * G, int ncount, int ecount, int *elist,
        int forward)
{
    int i, rval = 0;
    int *p;
    COLORadjnode *nodelist;

    COLORadjgraph_init (G);
    G->ncount = ncount;
    G->ecount = ecount;

    G->nodelist = COLOR_SAFE_MALLOC (G->ncount, COLORadjnode);
    COLORcheck_NULL (G->nodelist, "out of memory for G->nodelist");
    nodelist = G->nodelist;

    if (G->ecount) {
        G->adjspace = COLOR_SAFE_MALLOC (G->ecount, int);
        COLORcheck_NULL (G->adjspace, "out of memory for G->adjspace");
    }

    for (i = 0; i < ncount; i++) {
        nodelist[i].degree = 0;
    }

    for (i = 0; i < ecount; i++) {
        if (forward) {
            nodelist[elist[2*i]].degree++;
        } else {
            nodelist[elist[2*i+1]].degree++;
        }
    }

    p = G->adjspace;
    for (i = 0; i < ncount; i++) {
        nodelist[i].adj = p;
        p += nodelist[i].degree;
        nodelist[i].degree = 0;
    }

    for (i = 0; i < ecount; i++) {
        if (forward) {
            nodelist[elist[2*i]].adj[nodelist[elist[2*i]].degree++] = 
                                                               elist[2*i+1];
        } else {
            nodelist[elist[2*i+1]].adj[nodelist[elist[2*i+1]].degree++] = 
                                                               elist[2*i];
        }
    }

    /* Sort adj lists in backwards order */

    for (i = 0; i < ncount; i++) {
        if (nodelist[i].degree) {
            COLORutil_quicksort_reverse (nodelist[i].adj, nodelist[i].degree);
        }
    }

CLEANUP:

    if (rval) COLORadjgraph_free (G);
    return rval;
}
