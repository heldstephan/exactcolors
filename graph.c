#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "color.h"
#include "graph.h"

int COLORadjgraph_build (graph* G,
                         int ncount,int ecount, const int elist[])
{
    int rval = 0;
    int i;
    int *p;
    node *nodelist;

    G->nodelist = (node *) NULL;
    G->adjspace = (int *) NULL;
    G->ncount = ncount;
    G->ecount = ecount;

    G->nodelist = (node *) malloc (G->ncount * sizeof (node));
    if (!G->nodelist) {
        fprintf (stderr, "out of memory in build_graph\n");
        rval = 1; goto CLEANUP;
    }
    nodelist = G->nodelist;

    if (G->ecount) {
        G->adjspace = (int *) malloc (2 * G->ecount * sizeof (int));
        if (!G->adjspace) {
            fprintf (stderr, "out of memory in build_graph\n");
            rval = 1; goto CLEANUP;
        }
    }

    for (i = 0; i < ncount; i++) {
        nodelist[i].degree = 0;
    }

    for (i = 0; i < ecount; i++) {
        nodelist[elist[2*i]].degree++;
        nodelist[elist[2*i+1]].degree++;
    }

    p = G->adjspace;
    for (i = 0; i < ncount; i++) {
        nodelist[i].adj = p;
        p += nodelist[i].degree;
        nodelist[i].degree = 0;
    }

    for (i = 0; i < ecount; i++) {
        nodelist[elist[2*i]].adj[nodelist[elist[2*i]].degree++] =
                                                         elist[2*i+1];
        nodelist[elist[2*i+1]].adj[nodelist[elist[2*i+1]].degree++] =
                                                         elist[2*i];
    }

CLEANUP:

    if (rval) COLORadjgraph_free (G);
    return rval;
}

 void COLORadjgraph_free (graph *G)
{
    if (G->nodelist) free (G->nodelist);
    if (G->adjspace) free (G->adjspace);
}

static int comp_node_ids(const void* v1, const void* v2)
{
   int id1 = * (const int*) v1;
   int id2 = * (const int*) v2;

   return id1 - id2;
}

static void swap_nodes(int* v1, int* v2)
{
   int tmp = *v1;
   *v1     = *v2;
   *v2     = tmp;
}

static int unify_adjlist(int* adjlist,int degree, int* tmp_adjlist)
{
   int j;
   int new_degree = 0;

   if (degree) {
      tmp_adjlist[0] = adjlist[0];
      new_degree++;
      for (j = 1; j < degree; ++j) {
         if (adjlist[j] != adjlist[j-1]) {
            tmp_adjlist[new_degree++] = adjlist[j];
         }
      }
      for (j = 0; j < new_degree; ++j) {
            adjlist[j] = tmp_adjlist[j] ;
      }
   }
   return new_degree;
}

int COLORadjgraph_simplify(graph* G)
{
   int i,j;
   int rval = 0;
   int* tmp_adjlist = (int* ) NULL;

   assert(G);

   /* Create a sufficiently large working array.*/
   tmp_adjlist = (int*) malloc(G->ecount * sizeof(int));
   COLORcheck_NULL(tmp_adjlist,"Failed allocating tmp_adjlist");

   for (i = 0; i < G->ncount;++i) {
      int new_degree;
      int nloops = 0;

      qsort(G->nodelist[i].adj,G->nodelist[i].degree,sizeof(int),comp_node_ids);
      new_degree = unify_adjlist(G->nodelist[i].adj,G->nodelist[i].degree,
                                 tmp_adjlist);

      if(new_degree != G->nodelist[i].degree) {
         printf("Removed %d edge(s) from node %d.\n", 
                G->nodelist[i].degree - new_degree, i);
      }
      G->nodelist[i].degree = new_degree;
      
      for (j = 0; j < G->nodelist[i].degree; ++j) {
         if (G->nodelist[i].adj[j] == i) {
            nloops++;
            swap_nodes( & (G->nodelist[i].adj[j]), 
                        & (G->nodelist[i].adj[G->nodelist[i].degree - 1]) ); 
            --G->nodelist[i].degree;
            --j;
         }
      }
      if (nloops) {
         printf("Removed %d loop(s) from node %d.\n", nloops,i);
      }
            
   }
 CLEANUP:
   if (tmp_adjlist) {free(tmp_adjlist);}
   return rval;
}

int COLORadjgraph_extract_edgelist(int* ecount, int* elist[], const graph* G)
{
   int rval = 0;
   int i;
   *ecount = 0;
   if (*elist) {free(*elist);}

   for (i = 0; i < G->ncount;++i) {
      *ecount += G->nodelist[i].degree;
   }
   assert(*ecount % 2 == 0);
   /* elist of of size 2 * number of edges (== current *ecount).*/
   (*elist) = (int*) malloc( (*ecount) * sizeof(int));
   *ecount = 0;
   for (i = 0; i < G->ncount;++i) {
      int j;
      for (j = 0; j < G->nodelist[i].degree; ++j) {
         if (G->nodelist[i].adj[j] > i) {
            (*elist)[(*ecount) * 2]     = i;
            (*elist)[(*ecount) * 2 + 1] = G->nodelist[i].adj[j];
            (*ecount)++;
         }
      }
   }

   return rval;
}

void COLORadjgraph_sort_adjlists_by_id(graph* G)
{
   int i;
   for (i = 0; i < G->ncount;++i) {
      qsort(G->nodelist[i].adj,G->nodelist[i].degree,sizeof(int),comp_node_ids);
   }
}

