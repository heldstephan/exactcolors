#include <stdlib.h>
#include <stdio.h>
#include "graph.h"

#include "color.h"


static void color_node (graph *G, int n);
static void perm_quicksort (int *perm, int *len, int n);

int COLORgreedy (int ncount, int ecount, int *elist, int *ncolors,
        COLORset **colorclasses)
{
    int rval = 0;
    int *degree = (int *) NULL;
    int *perm = (int *) NULL;
    int i, k, c;
    graph G;
    COLORset *csets = (COLORset *) NULL;

    printf ("COLORgreedy(%d,%d) ...\n", ncount, ecount);
    fflush (stdout);

    *ncolors = 0;
    *colorclasses = (COLORset *) NULL;
    
    rval = COLORadjgraph_build (&G, ncount, ecount, elist);
    if (rval) {
        fprintf (stderr, "build_graph failed\n");
        goto CLEANUP;
    }

    degree = (int *) malloc (ncount * sizeof (int));
    if (!degree) {
        fprintf (stderr, "out of memory for degree\n");
        rval = 1;  goto CLEANUP;
    }
    perm = (int *) malloc (ncount * sizeof (int));
    if (!perm) {
        fprintf (stderr, "out of memory for perm\n");
        rval = 1;  goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) {
        degree[i] = G.nodelist[i].degree;
        perm[i] = i;
    }

    perm_quicksort (perm, degree, ncount);
    for (i = 0; i < ncount; i++) G.nodelist[i].color = -1;

    for (i = 0; i < ncount; i++) {
        color_node (&G, perm[i]);
    } 

    k = 0;
    for (i = 0; i < ncount; i++) {
        if (G.nodelist[i].color > k) {
            k = G.nodelist[i].color;
        }
    }
    k++;
    *ncolors = k;
    csets = (COLORset *) malloc (k * sizeof (COLORset));
    if (!csets) {
        fprintf (stderr, "out of memory for csets\n");
        rval = 1;  goto CLEANUP;
    }

    for (i = 0; i < k; i++) {
        COLORinit_set (&csets[i]);
    }
    for (i = 0; i < ncount; i++) {
        csets[G.nodelist[i].color].count++;
    }
    for (i = 0; i < k; i++) {
        csets[i].members = (int *) malloc (csets[i].count * sizeof (int));
        if (!csets[i].members) {
            fprintf (stderr, "out of memory for csets members\n");
            rval = 1; goto CLEANUP;
        }
        csets[i].count = 0;
    }

    for (i = 0; i < ncount; i++) {
        c = G.nodelist[i].color;
        csets[c].members[csets[c].count] = i;
        csets[c].count++;
    }
    *colorclasses = csets;

CLEANUP:

    if (rval) {
        if (csets) {
            for (i = 0; i < k; i++) COLORfree_set (&csets[i]);
            free (csets);
        }
    }
    if (degree) free (degree);
    if (perm) free (perm);
    COLORadjgraph_free (&G);
    return rval;
}

static void color_node (graph *G, int n)
{
    int i, color = -1;
    node *p = &G->nodelist[n];

    for (i = 0; i < p->degree; i++) {
        if (G->nodelist[p->adj[i]].color > color) {
            color = G->nodelist[p->adj[i]].color;
        }
    }
    p->color = color+1;
}

static void perm_quicksort (int *perm, int *len, int n)
{
    int i, j, temp, t;

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

    perm_quicksort (perm, len, j);
    perm_quicksort (perm + i, len, n - i);
}

void COLORinit_set (COLORset *s) 
{
    if (s) {
        s->members = (int *) NULL;
        s->count = 0;
    }
}

void COLORfree_set (COLORset *s) 
{
    if (s && s->members) {
        free (s->members);
        s->members = (int *) NULL;
        s->count = 0;
    }
}

void COLORfree_sets(COLORset** sets,int* nsets)
{
   int i;
   if (*sets) {
      for (i = 0; i < *nsets; i++) COLORfree_set (& (*sets)[i]);
      free (*sets);
   }
   *sets  = (COLORset*) NULL;
   *nsets = 0;
}
