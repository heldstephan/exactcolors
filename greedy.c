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
#include <assert.h>
#include "graph.h"

#include "color.h"


static void color_node (COLORadjgraph *G, int n);
static void update_satdegree_for_neighbors (COLORadjgraph *G, int colored_node,
                                            int freedegree[], int satdegree[], int tmp_colors[]);


int COLORgreedy (int ncount, int ecount, int *elist, int *ncolors,
        COLORset **colorclasses)
{
    int rval = 0;
    int *degree = (int *) NULL;
    int *perm = (int *) NULL;
    int i, k, c;
    COLORadjgraph G;
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

    COLORutil_perm_quicksort (perm, degree, ncount);
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

int COLORdsatur (int ncount, int ecount, int *elist, int *ncolors,
                 COLORset **colorclasses)
{
   int rval = 0;
   int *freedegree = (int *) NULL;
   int *satdegree = (int *) NULL; /* saturation degree.*/

   int *tmpcolors = (int *) NULL;

   int i, k, c;
   COLORadjgraph G;
   COLORset *csets = (COLORset *) NULL;

   printf ("COLORgreedy(%d,%d) ...\n", ncount, ecount);
   fflush (stdout);

   rval = COLORadjgraph_build (&G, ncount, ecount, elist);
   if (rval) {
      fprintf (stderr, "build_graph failed\n");
      goto CLEANUP;
   }

   freedegree = (int *) malloc (ncount * sizeof (int));
   if (!freedegree) {
      fprintf (stderr, "out of memory for freedegree\n");
      rval = 1;  goto CLEANUP;
   }

   satdegree = (int *) malloc (ncount * sizeof (int));
   if (!satdegree) {
      fprintf (stderr, "out of memory for satdegree\n");
      rval = 1;  goto CLEANUP;
   }

   tmpcolors = (int *) malloc (ncount * sizeof (int));
   if (!tmpcolors) {
      fprintf (stderr, "out of memory for tmpcolors\n");
      rval = 1;  goto CLEANUP;
   }

   for (i = 0; i < ncount; i++) {
      freedegree[i] = G.nodelist[i].degree;
      satdegree[i] = 0;
      G.nodelist[i].color = -1;
   }

   for (i = 0; i < *ncolors; ++i) {
      int j;
      for (j = 0; j < (*colorclasses)[i].count; ++j) {
         int v = (*colorclasses)[i].members[j];
         color_node (&G, v);
         update_satdegree_for_neighbors (&G, v,
                                         freedegree, satdegree, tmpcolors);
      }
   }

   COLORfree_sets(colorclasses,ncolors);

   for (i = 0; i < ncount - *ncolors; i++) {
      int maxsatdegree      = -1;
      int maxfreedegree     = -1;
      int maxsatdegree_node = -1;
      int j;

      for (j = 0; j < ncount; j++) {

         if (G.nodelist[j].color == -1) {
            if ( (satdegree[j] > maxsatdegree) ||
                 ( (satdegree[j] ==  maxsatdegree) &&
                   (freedegree[j] > maxfreedegree) ) )
            {
               maxsatdegree       = satdegree[j];
               maxfreedegree      = freedegree[j];
               maxsatdegree_node  = j;
            }
         }
      }
      if (maxsatdegree_node > -1) {
         color_node (&G, maxsatdegree_node);
         update_satdegree_for_neighbors (&G, maxsatdegree_node,
                                         freedegree, satdegree, tmpcolors);
      }
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
   COLOR_IFFREE(freedegree,int);
   COLOR_IFFREE(satdegree,int);
   COLOR_IFFREE(tmpcolors,int);

   COLORadjgraph_free (&G);
   return rval;
}

static void update_satdegree_for_neighbors (COLORadjgraph *G, int colored_node,
                                            int freedegree[], int satdegree[], int tmpcolors[])
{
   int j;
   for (j = 0;j < G->nodelist[colored_node].degree; ++j) {
      int v = G->nodelist[colored_node].adj[j];
      int k;
      freedegree[v]--;
      satdegree[v] = 0;
      for (k = 0; k < G->ncount;++k) {
         tmpcolors[k] = 0;
      }
      for (k = 0; k < G->nodelist[v].degree;++k) {
         int w = G->nodelist[v].adj[k];
         int wcolor = G->nodelist[w].color;
         if (wcolor>-1 && tmpcolors[wcolor] == 0) {
            tmpcolors[wcolor] = 1;
            satdegree[v]++;
         }
      }
   }
 }


static void color_node (COLORadjgraph *G, int n)
{
    int i, color = 0;
    COLORadjnode *p = &G->nodelist[n];
    int failed;
    do {
     failed = 0;
       for (i = 0; !failed && i < p->degree; i++) {
          if (G->nodelist[p->adj[i]].color == color) {
             ++color;
             failed = 1;
          }
       }
    } while (failed);
    p->color = color;
}

int COLORtransform_into_maximal (int ncount, int ecount, int *elist, int ncolors,
                                 COLORset *cclasses)
{
   int rval = 0;
   int *color_incidence= (int *) NULL;
   COLORadjgraph G;
   int i,c;


   rval = COLORadjgraph_build (&G, ncount, ecount, elist);
   COLORcheck_rval(rval, "Failed in COLORadjgraph_build.");

   color_incidence = COLOR_SAFE_MALLOC(ncount,int);
   COLORcheck_NULL(color_incidence,"Failed to allocate color_incidence.");

   for (c = 0; c < ncolors; ++c) {
      int new_count = cclasses[c].count;
      for (i = 0; i < ncount; ++i) color_incidence[i] = 0;
      for (i = 0; i < cclasses[c].count; ++i) {
         color_incidence[cclasses[c].members[i]] = 1;
      }

      for (i = 0; i < ncount; ++i) {
         int j;
         int free_node = 1;
         if (color_incidence[i]) continue;
         for (j = 0 ; j < G.nodelist[i].degree;++j) {
            if (color_incidence[G.nodelist[i].adj[j]]) {
               free_node = 0;
               j = G.nodelist[i].degree;
            }
         }
         if (free_node) {
            color_incidence[i] = 1;
            new_count++;
         }
      }

      if (new_count > cclasses[c].count) {
         COLOR_FREE(cclasses[c].members,int);
         cclasses[c].members = COLOR_SAFE_MALLOC(new_count,int);
         COLORcheck_NULL(cclasses[c].members, "Failed to allocate cclasses[c].members");
         cclasses[c].count = new_count;
         new_count = 0;
         for (i = 0; i < ncount; ++i) {
            if (color_incidence[i]) {
               cclasses[c].members[new_count++] = i;
            }
         }
         assert(new_count = cclasses[c].count);
      }
   }
 CLEANUP:
   return rval;
}


int COLORtransform_into_coloring(int ncount, int *ncolors,
                                 COLORset **colorclasses)
{
   int rval = 0;
   int* node_color   = (int*) NULL;

   int nnewcolors = 0;
   COLORset* newcolors = (COLORset*) NULL;
   int i,c;

   node_color = COLOR_SAFE_MALLOC(ncount,int);
   COLORcheck_NULL(node_color, "Failed to allocate node_color");


   for (i = 0; i< ncount; ++i) {
      node_color[i]   = -1;
   }
   for (c = 0; c < *ncolors; ++c) {
      int nc = 0;
      for (i = 0; i < (*colorclasses)[c].count; ++i) {
         if (node_color[ (*colorclasses)[c].members[i]] == -1) {
            node_color[(*colorclasses)[c].members[i]] = nnewcolors;
            nc++;
         } else {
            // current node is already covered.
         }
      }
      if (nc) {
         nnewcolors++;
      }
   }

   COLORfree_sets (colorclasses, ncolors);

   newcolors = COLOR_SAFE_MALLOC(nnewcolors, COLORset);
   COLORcheck_NULL(newcolors, "Failed to allocate newcolors");

   for (i = 0; i < nnewcolors; i++) {
      COLORinit_set (newcolors + i);
   }
   for (i = 0; i < ncount; i++) {
      newcolors[node_color[i]].count++;
   }

   for (c = 0; c < nnewcolors; c++) {
      newcolors[c].members = COLOR_SAFE_MALLOC(newcolors[c].count,int);
      COLORcheck_NULL(newcolors[c].members,"Failed to allocate newcolors[c].members.");

      newcolors[c].count = 0;
   }

   for (i = 0; i < ncount; i++) {
      int c = node_color[i];
      newcolors[c].members[newcolors[c].count] = i;
      newcolors[c].count++;
   }

   *colorclasses = newcolors;
   *ncolors      = nnewcolors;
   newcolors  = (COLORset*) NULL;
   nnewcolors = 0;


CLEANUP:

   COLORfree_sets (&newcolors, &nnewcolors);
   COLOR_IFFREE(node_color, int);

   return rval;
}
