#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "color.h"
#include "graph.h"
#include "color_defs.h"

int COLORadjgraph_build (graph* G, int ncount, int ecount, const int elist[])
{
    int rval = 0;
    int i;
    int *p;
    node *nodelist;

    COLORadjgraph_init (G);
    G->ncount = ncount;
    G->ecount = ecount;

    G->nodelist = COLOR_SAFE_MALLOC (G->ncount, node);
    COLORcheck_NULL (G->nodelist, "out of memory for G->nodelist");
    nodelist = G->nodelist;

    if (G->ecount) {
        G->adjspace = COLOR_SAFE_MALLOC (2 * G->ecount, int);
        COLORcheck_NULL (G->adjspace, "out of memory for G->adjspace");
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


int  COLORadjgraph_copy(graph* Gdst, const graph* Gsrc)
{
   /* This is a fast too implement version using ecxisting functions.
      Copying graphs couls be done faster.
      S. Held.
   */
   int rval = 0;
   int* elist = (int*) NULL;
   int ecount = 0;
   
   rval = COLORadjgraph_extract_edgelist(&ecount,&elist,Gsrc);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_extract_edgelist");
   
   rval = COLORadjgraph_build(Gdst,Gsrc->ncount,ecount,elist);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_build");
   
 CLEANUP:
   if (rval) COLORadjgraph_free(Gdst);
   if (elist) free(elist);
   return rval;   
}

int  COLORadjgraph_build_complement(graph* Gc, const graph* G)
{
   int rval = 0;
   int v_i, a_i,na;
   int*  elist = (int*) NULL;
   int  ecount = 0, ecount_chk = 0;
   
   rval = COLORadjgraph_copy(Gc, G);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_copy.");

   /* Simplify will also sort the adjaceny lists.*/
   rval = COLORadjgraph_simplify(Gc);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_simplify.");
   
   ecount_chk = (Gc->ncount * (Gc->ncount - 1) )/ 2 - Gc->ecount;
   

   elist = COLOR_SAFE_MALLOC (2*ecount_chk, int);
   COLORcheck_NULL(elist,"Failed to allocate elist.");
   
   
   ecount = 0;  
   for (v_i = 0; v_i < Gc->ncount; ++ v_i) {
      node* v = &(Gc->nodelist[v_i]);
      int a = -1;
      a_i  = 0;
      for (na = v_i + 1; na < Gc->ncount; ++ na) {
         while (a_i < v->degree && a < na) {
            a = v->adj[a_i];
            ++a_i;
         }
         if (na != a) {
            elist[2*ecount]   = v_i;
            elist[2*ecount+1] = na;
            ++ecount;
         }
      }
   }
   assert(ecount == ecount_chk);

   COLORadjgraph_free(Gc);
   
   rval = COLORadjgraph_build(Gc, G->ncount,ecount,elist);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_build");
 CLEANUP:
   if (rval) {
      COLORadjgraph_free(Gc);
   }
   if (elist) free(elist);
   return rval;
}


void COLORadjgraph_init (graph *G)
{
   if (G) {
      G->nodelist = (node *) NULL;
      G->adjspace = (int *) NULL;
      G->ncount = 0;
      G->ecount = 0;
   }
}

void COLORadjgraph_free (graph *G)
{
    if (G->nodelist) free (G->nodelist);
    if (G->adjspace) free (G->adjspace);
        
    G->nodelist = (node*) NULL;
    G->adjspace = (int*)  NULL;
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
   int ncount, ecount;

   assert(G);

   /* Create a sufficiently large working array.*/
   tmp_adjlist = COLOR_SAFE_MALLOC (G->ecount, int);
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

   /* Re-allocate graph to generate correct ecount.*/
   if (tmp_adjlist) {
      free(tmp_adjlist);
      tmp_adjlist = (int*) NULL;
   }
   ncount = G->ncount;
   rval = COLORadjgraph_extract_edgelist(&ecount,&tmp_adjlist,G);
   COLORcheck_rval(rval, "Failed in COLORadjgraph_extract_edgelist");

   COLORadjgraph_free(G);
   
   rval = COLORadjgraph_build(G, ncount,ecount,tmp_adjlist);
   COLORcheck_rval(rval, "Failed in COLORadjgraph_build");
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
   (*elist) = COLOR_SAFE_MALLOC ((*ecount), int);
   COLORcheck_NULL (*elist, "out of memory for elist");
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

CLEANUP:

   return rval;
}

void COLORadjgraph_sort_adjlists_by_id(graph* G)
{
   int i;
   for (i = 0; i < G->ncount;++i) {
      qsort(G->nodelist[i].adj,G->nodelist[i].degree,sizeof(int),comp_node_ids);
   }
}

int  COLORadjgraph_delete_unweighted(graph* G, 
                                     int** new_nweights,
                                     const int nweights[])
{
   int rval = 0;
   int* nmap = (int*) NULL;
   int* newelist = (int*) NULL;
   int i,a_i;
   int ncount = 0;
   int ecount = 0;

   nmap = COLOR_SAFE_MALLOC (G->ncount, int);
   COLORcheck_NULL(nmap,"Failed to allocate nmap");

   newelist = COLOR_SAFE_MALLOC (2*G->ecount, int);
   COLORcheck_NULL(nmap,"Failed to allocate newelist");

   for (i = 0; i < G->ncount; ++i) {
      if (nweights[i] == 0) {
         nmap[i] = -1;
         G->nodelist[i].degree = 0;
      } else {
         int a = 0;
         nmap[i] = ncount++;
         
         for (a_i = 0; a_i < G->nodelist[i].degree && a < i; ++a_i) {
            a = G->nodelist[i].adj[a_i];
            if (a < i && nmap[a] != -1) {
               newelist[2*ecount] = nmap[a];
               newelist[2*ecount+1] = nmap[i];
               ++ ecount;
            }
         }
      }
   }
   *new_nweights = COLOR_SAFE_MALLOC (ncount, int);
   COLORcheck_NULL(*new_nweights,"Failed to allocate nmap");

   for (i = 0; i < G->ncount; ++i) {
      int ni = nmap[i];
      if (ni != -1) {
         (*new_nweights)[ni] = nweights[i];
      }
   }
              
   COLORadjgraph_free(G);

   rval = COLORadjgraph_build(G,ncount,ecount,newelist);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_build");

   
   printf("Reduced graph has %d nodes and %d edges.\n",
          ncount,ecount);

   

 CLEANUP:
   if (nmap) free(nmap);
   if (newelist) free(newelist);
   if (rval) {
      if (*new_nweights) free(*new_nweights);
      *new_nweights = (int*) NULL;
   }
   return rval;
}

int COLORread_dimacs (char *f, int *pncount, int *pecount, int **pelist,
       int **pnweights)
{
    int rval = 0;
    int ncount, ecount, icount = 0, haveprob = 0;
    int i, end0, end1, n, len;
    int *elist = (int *) NULL;
    int *nweights = (int *) NULL;
    char buf[256], *p;
    FILE *in = (FILE *) NULL;
    graph G; /* used to simplify graph.*/

    G.nodelist = (node *) NULL;
    G.adjspace = (int *) NULL;
    G.ncount = 0;
    G.ecount = 0;

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
            if ( strcmp(data,"edge") && strcmp(data,"edges") &&
                                        strcmp(data,"col") ) {
                fprintf (stderr, "ERROR in Dimacs file -- not an edge file\n");
                rval = 1;  goto CLEANUP;
            }
            data = strtok(NULL,delim);
            sscanf (data, "%d", &ncount);
            data = strtok(NULL,delim);
            sscanf (data, "%d", &ecount);

            printf ("Number of Nodes: %d\n", ncount);
            printf ("Number of Edges: %d\n", ecount);
            elist = COLOR_SAFE_MALLOC (2*ecount, int);
            COLORcheck_NULL (elist, "out of memory for elist");
            nweights = COLOR_SAFE_MALLOC (ncount, int);
            COLORcheck_NULL (nweights, "out of memory for nweights");
            for (i = 0; i < ncount; i++) nweights[i] = 0;
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
            nweights[n-1] = len;
        }
    }

    rval = COLORadjgraph_build(&G, ncount,icount,elist);
    COLORcheck_rval(rval,"COLORadjgraph_build failed");                                     
    rval = COLORadjgraph_simplify(&G);
    COLORcheck_rval(rval,"COLORadjgraph_simplify failed");                                     
    COLORadjgraph_extract_edgelist(&icount, &elist,&G);
    COLORcheck_rval(rval,"COLORadjgraph_extract_edgelist");                                     
    *pncount = ncount;
    /* Some col-instances are buggy => reduce # edges to icount*/
    *pecount = icount; 
    *pelist = elist;
    if (pnweights) {
        *pnweights = nweights;
    } else {
        COLOR_IFFREE (nweights, int);
    }

CLEANUP:

    COLORadjgraph_free(&G);
    if (rval) {
        COLOR_IFFREE (elist, int);
        COLOR_IFFREE (nweights, int);
    }
    if (in) fclose (in);

    return rval;
}

int COLORedge_stat(const graph* G)
{
   int rval = 0;
   int i;
   int* degreecnt = (int*) NULL;
   
   degreecnt = COLOR_SAFE_MALLOC (G->ncount, int);
   COLORcheck_NULL(degreecnt,"Failed to allocate degreecnt");

   for (i = 0; i < G->ncount; ++i) {
      degreecnt[i] = 0;
   }

   for (i = 0; i < G->ncount; ++i) {
      ++(degreecnt[G->nodelist[i].degree]);
   }

   for (i = 0; i < G->ncount; ++i) {
      if (degreecnt[i]) {
         printf("DEG %d NUM %d\n",i, degreecnt[i]);
      }
   }
   
 CLEANUP:
   if (degreecnt) free(degreecnt);
   return rval;
}

int  COLORgraph_print(int ecount, const int elist[])
{
   int i;
   for (i = 0; i < ecount; ++i) {
      printf("e %d %d\n",elist[2*i], elist[2*i+1]);
   }
   return 0;
}
