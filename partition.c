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
#include <string.h>
#include <limits.h>
#include <getopt.h>

#include "graph.h"
#include "heap.h"
#include "color.h"


char* graphfile = (char*) NULL;
int   partitionsize = INT_MAX;

static int parseargs (int ac, char **av);
static void usage (char *f);
static int get_problem_name(char* pname,const char* efname);

static int COLORadjgraph_isolate_vertex(COLORadjgraph* G, int v_i,
				 COLORNWTHeap* heap, int href[])
{
   int rval = 0;
   int j;
   COLORadjnode* v = G->nodelist + v_i;

   for (j = 0; j < v->degree; ++j) {
      COLORadjnode* w = G->nodelist + v->adj[j];
      int k;
      if (heap) {
	 COLORNWTheap_decrease_key(heap,href[v->adj[j]],w->degree -1);
      }
      for (k = 0; k < w->degree; ++k) {
	 if (w->adj[k] == v_i) {
	    printf("c Deleting %d %d\n", v_i + 1, v->adj[j] + 1);
	    if (k + 1 < w->degree) {
	       memmove(w->adj + k , w->adj + k + 1,
		       (w->degree - k -1) * sizeof(int));
	    }
	    w->degree--;
	    k = w->degree;
	 }
      }
   }

   v->degree = 0;

   return rval;
}

COLOR_MAYBE_UNUSED
static
int COLORadjgraph_incr_deled_smallest_degree_vertex(COLORadjgraph* Gd, 
						    const COLORadjgraph* Gs, 
						    int lpartitionsize)
{
   int           rval = 0;
   int*          heapref = (int*) NULL;
   int*          nodeindices = (int*) NULL;
   int           ndeletednodes = 0;
   int           ndeletededges = 0;
   COLORNWTHeap* heap = (COLORNWTHeap*) NULL;
   int i;
   int           ndelete = Gs->ncount - lpartitionsize;
   int* i_ptr;

   rval =   COLORadjgraph_copy(Gd, Gs);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_copy");

   /* rval = COLORadjgraph_extract_edgelist(&ecount, &elist, Gd); */
   /* COLORcheck_rval(rval,"Failed in COLORadjgraph_extract_edgelist"); */

   /* rval =   COLORgraph_print(ecount, elist); */
   /* COLORcheck_rval(rval,"Failed in COLORgraph_print"); */

   rval = COLORNWTheap_init(&heap, Gs->ncount);
   COLORcheck_rval(rval,"Failed in COLORNWTheap_init");

   heapref = COLOR_SAFE_MALLOC(Gd->ncount, int);
   COLORcheck_NULL(heapref,"Failed to allocate heapref");

   nodeindices = COLOR_SAFE_MALLOC(Gd->ncount, int);
   COLORcheck_NULL(heapref,"Failed to allocate nodeindices");


   for (i = 0; i < Gd->ncount; ++i) {
      COLORadjnode* node = Gd->nodelist + i;
      COLORNWT      key = (COLORNWT) node->degree;
      int*          pos = &(heapref[i]);
      nodeindices[i] = i;
      COLORNWTheap_insert(heap, pos,key,nodeindices + i);
   }
   
   while ( (ndeletednodes < ndelete) &&(i_ptr = (int*) COLORNWTheap_min(heap) ) ) {
      i = *i_ptr;
      ndeletednodes++;
      ndeletededges += Gd->nodelist[i].degree;
      rval = COLORadjgraph_isolate_vertex(Gd,i,heap,heapref);
      COLORcheck_rval(rval, "Failed in COLORadjgraph_isolate_vertex");
      
      /* rval = COLORadjgraph_extract_edgelist(&ecount, &elist, Gd); */
      /* COLORcheck_rval(rval,"Failed in COLORadjgraph_extract_edgelist"); */
      /* printf("edgelist after isolating %d.\n",i); */

      /* rval =   COLORgraph_print(ecount, elist); */
      /* COLORcheck_rval(rval,"Failed in COLORgraph_print"); */

   }

 CLEANUP:
   COLORNWTheap_free(heap);
   return rval;
}

COLOR_MAYBE_UNUSED
static
int COLORadjgraph_greedy_partition_1(COLORadjgraph* Gd, 
				     const COLORadjgraph* Gs, 
				     int lpartitionsize)
{
   int rval = 0;
   int* perm = (int*) NULL;
   int* iperm = (int*) NULL;
   int* len  = (int*) NULL;
   int i;
   int ncount     =  Gs->ncount;
   int nbest_half = (lpartitionsize + 1) / 2;
   int nremaining =  ncount - nbest_half;

   rval =   COLORadjgraph_copy(Gd, Gs);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_copy");

   perm = COLOR_SAFE_MALLOC(ncount,int);
   COLORcheck_NULL(perm,"Failed to allocate perm.\n");

   iperm = COLOR_SAFE_MALLOC(ncount,int);
   COLORcheck_NULL(iperm,"Failed to allocate iperm.\n");

   len = COLOR_SAFE_MALLOC(ncount,int);
   COLORcheck_NULL(len,"Failed to allocate len.\n");
   
   for (i=0; i < ncount; ++i) {
      perm[i] = i;
      len[i] = Gd->nodelist[i].degree;
   }

   COLORutil_perm_quicksort (perm,len,ncount);

   for (i=0; i < ncount; ++i) {
      iperm[perm[i]] = i;
   }
   for (i=0; i < nremaining; ++i) {
      COLORadjnode* w = Gd->nodelist + perm[i];
      int j;
      len[perm[i]] = 0;
      for (j = 0; j < w->degree; ++j) {
	 int k = w->adj[j];
	 if (iperm[k] > nremaining) {
	    len[perm[i]]++;
	 }
      }
   }
   COLORutil_perm_quicksort (perm,len,nremaining);
      
   for (i=0; i < nremaining; ++i) {
      iperm[perm[i]] = i;
   }
   for( i = 0; i < ncount - lpartitionsize; i++) {
      rval = COLORadjgraph_isolate_vertex(Gd, iperm[i],NULL,NULL);
      COLORcheck_rval(rval,"Failed in COLORadjgraph_isolate_vertex");
   }


 CLEANUP:
   COLOR_IFFREE(perm,int);
   return rval;
}
static int parseargs (int ac, char **av)
{
    int c;
    int rval = 0;

    while ((c = getopt (ac, av, "dp:")) != EOF) {
        switch (c) {
        case 'd':
	   COLORset_dbg_lvl(COLORdbg_lvl() + 1);
            break;
        case 'p':
	   partitionsize = atoi(optarg);
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
    fprintf (stderr, "   -d    turn on debugging\n");
    fprintf (stderr, "   -p n  targeted max. partition size.\n");
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

   return 0;
}




int main (int ac, char **av)
{
    char pname[256] = "";
    int rval = 0;
    int ncount, ecount;
    int *elist = (int *) NULL;
    int *wlen = (int *) NULL;
    int *new_ids = (int *) NULL;
    
    int ndelete;
    COLORadjgraph G1, G2;
    int i, new_ncount;

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    get_problem_name (pname, graphfile);

    printf("c Graph generated from %s with <= %d non-isolated vertices.\n",
	   pname, partitionsize);

    rval = COLORread_dimacs (graphfile, &ncount, &ecount, &elist, &wlen);
    COLORcheck_rval (rval, "COLORread_dimacs failed");

    rval = COLORadjgraph_build (&G1, ncount, ecount, elist);
    COLORcheck_rval (rval, "COLORadjgraph_build failed");

    new_ids = COLOR_SAFE_MALLOC(ncount, int);
    COLORcheck_NULL(new_ids, "Failed to allocate new_ids");
    for (i = 0; i < ncount; ++i) { new_ids[i] = i;}

    new_ncount = ncount;
    ndelete = ncount - partitionsize;
    if (ndelete > 0) {

/*        rval = COLORadjgraph_greedy_partition_1 (&G2, &G1,partitionsize); */
/*        COLORcheck_rval (rval, "COLORadjgraph_build failed"); */

       COLORadjgraph_incr_deled_smallest_degree_vertex(&G2, 
                                                       &G1,
                                                       partitionsize);


       rval = COLORadjgraph_extract_edgelist(&ecount, &elist, &G2);
       COLORcheck_rval(rval,"Failed in COLORadjgraph_extract_edgelist");
    

       new_ncount = 0;
       for (i = 0; i < ncount; ++i) { 
          COLORadjnode* v = G2.nodelist + i;
          if (v->degree){
             new_ids[i] = new_ncount++;
             printf("c node %d renamed to %d\n", i + 1, new_ids[i] + 1);
          } else {
             new_ids[i] = -1;
          }
       }
    }

    printf("p edge %d %d\n", new_ncount,ecount);
    for (i = 0; i < ecount; ++i) {
       printf("e %d %d\n",new_ids[elist[2*i]] + 1,new_ids[elist[2*i+1]] + 1);
    }
 CLEANUP:
    return rval;
}
