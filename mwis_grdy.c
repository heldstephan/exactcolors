#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <assert.h>

#include <valgrind/valgrind.h>

#include "lp.h"

#include "graph.h"
#include "heap.h"

#include "mwis.h"


typedef struct COLORclasses {
   int       cnt;
   int       allocnt;
   COLORset* sets;
} COLORclasses;


static void COLORclasses_init(COLORclasses* cclasses)
{
   cclasses->sets =  (COLORset*) NULL;
   cclasses->cnt     = 0;
   cclasses->allocnt = 0;
}



static void COLORclasses_reset_without_free(COLORclasses* cclasses)
{
   while(cclasses->cnt){
      cclasses->cnt--;
      cclasses->sets[cclasses->cnt].members = (int*) NULL;
      cclasses->sets[cclasses->cnt].count = 0;
   }
   assert(cclasses->cnt == 0);
}


static int COLORclasses_expand(COLORclasses* classes)
{
   int rval = 0;

   COLORset* oldsets = classes->sets;
   classes->allocnt = 1.5 * classes->allocnt + 1;

   classes->sets = (COLORset*) COLOR_SAFE_MALLOC(classes->allocnt,COLORset);
   COLORcheck_NULL(classes->sets,"Failed to allocatete classes->sets");

   memcpy(classes->sets,oldsets,classes->cnt * sizeof(COLORset));

 CLEANUP:
   if (oldsets) free(oldsets);
   if (rval) {
      if (classes->sets) free(classes->sets);
   }
   return rval;
}


/**
 The implementation for the local search data structures are inspired
 by corresponding data structures for the unweighted case from
 following paper:

 Diogo Viera Andrade, Mauricio G. C. Resende, Renato Fonseca
 F. Werneck: Fast Local Search for the Maximum Independent Set
 Problem.  WEA 2008, 220-234.
*/
typedef struct soldata_t soldata;
struct soldata_t {
   int*  nperm;   /* node permutation.*/
   int*  inperm;  /* inverse node permutation.*/
   int*  tightness;
   int   ncount;
   int   solcount;
   int   freecount;
   const graph*   G;
   const COLORNWT*  nweights;
   COLORNWT    cutoff;
   COLORclasses cclasses;

   /* Following arrays are used in diverse subroutines
      and kept in soldata avoid too many mallocs.
   */
   int*      sort_work;
   COLORNWT* grdy_nweights;
   int*      work_marker;
   int*      work_path;
   int*      nodestack;

   int*          heapref;
   COLORNWTHeap* heap;

   int (* greedy_improvement)(soldata*,int);

   COLORlp *lp;
   double*  dbl_nweights;
};

struct _MWISls_env {
   graph    G;
   soldata sol;
};


static void remove_from_soldata(soldata* sol, int v);
static void set_ptrs_to_zero(soldata* sol);


static void clean_soldata(soldata* sol) {
   if (sol->nperm)       free(sol->nperm);
   if (sol->inperm)      free(sol->inperm);
   if (sol->tightness)   free(sol->tightness);
   if (sol->sort_work)   free(sol->sort_work);
   if (sol->grdy_nweights)    free(sol->grdy_nweights);
   if (sol->work_marker) free(sol->work_marker);
   if (sol->work_path)   free(sol->work_path);
   if (sol->nodestack)   free(sol->nodestack);
   if (sol->heapref)     free(sol->heapref);
   if (sol->cclasses.sets) free(sol->cclasses.sets);
   if (sol->heap)
      COLORNWTheap_free(sol->heap);

   set_ptrs_to_zero(sol);
}


/* static int build_inverse_perm(soldata* sol) */
/* { */
/*    int i = 0; */
/*    for (i = 0; i < sol->freecount; ++i) { */
/*       int j = sol->solcount + i; */
/*       sol->inperm[sol->nperm[j]] = j; */
/*    } */
/*    return 0; */
/* } */

static int is_free(soldata* sol,int i)
{
   return (sol->solcount <= sol->inperm[i] &&
           sol->inperm[i] < sol->solcount + sol->freecount);
}

static int is_basis(soldata* sol,int i)
{
   return (sol->inperm[i] < sol->solcount);
}

static int swap_nodelist(soldata* sol,int i,int j)
{
   if (i == j) return 0;
   {
      int perm_i = sol->inperm[i];
      int perm_j = sol->inperm[j];

      int tmp_perm_i = sol->nperm[perm_i];      assert(tmp_perm_i == i);
      sol->nperm[perm_i] = sol->nperm[perm_j];
      sol->nperm[perm_j] = tmp_perm_i;

      sol->inperm[j] = perm_i;
      sol->inperm[i] = perm_j;
   }
   return 0;
}

static int tighten_neighbours(soldata* sol,int i)
{
   int j;
   node* nodelist = sol->G->nodelist;
   for( j = 0; j < nodelist[i].degree; ++j)
   {
      int k = nodelist[i].adj[j];
      if (COLORdbg_lvl() > 2) {
         printf("Increasing %d from %d to %d\n",
                k,sol->tightness[k],sol->tightness[k]+1);
      }

      sol->tightness[k]++;
      assert(!is_basis(sol,k));
      if (is_free(sol,k)) {
         int l = sol->nperm[sol->solcount + sol->freecount - 1];
         swap_nodelist(sol,k,l);
         sol->freecount--;
      }
   }
   return 0;
}



static int relax_neighbours(soldata* sol,int i)
{
   int j;
   node* nodelist = sol->G->nodelist;
   for( j = 0; j < nodelist[i].degree; ++j)
   {
      int k = nodelist[i].adj[j];
      if (COLORdbg_lvl() > 2) {
         printf("Decreasing %d from %d to %d\n",
                k,sol->tightness[k],sol->tightness[k]-1);
      }
      sol->tightness[k]--;
      if (sol->tightness[k] < 0) {
         printf("tightness if %d dropped to %d\n",
                k, sol->tightness[k]);
      }
      assert(sol->tightness[k] >= 0);
      /* k moves from non-free to free.*/
      if (sol->tightness[k] == 0 && !is_basis(sol,k)) {
         int l = sol->nperm[sol->solcount + sol->freecount];
         swap_nodelist(sol,k,l);
         sol->freecount++;
      }
   }
   return 0;
}

static int add_iff_free(soldata* sol,int i)
{
   if (is_free(sol,i)) {

      if (COLORdbg_lvl() > 1) {
         printf("Adding %d\n",i);
      }
      assert(sol->tightness[i] == 0);
      swap_nodelist(sol,i, sol->nperm[sol->solcount]);
      sol->solcount++;
      sol->freecount--;
      tighten_neighbours(sol,i);
      return 1;
   }
   return 0;
}

static void perm_nwt_rquicksort (int *perm, const COLORNWT *len, int n)
{
   int i, j, temp;
   COLORNWT t;

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

   perm_nwt_rquicksort (perm, len, j);
   perm_nwt_rquicksort (perm + i, len, n - i);
}

static void set_ptrs_to_zero(soldata* sol)
{
   sol->nperm         =      (int*) NULL;
   sol->inperm        =      (int*) NULL;
   sol->tightness     =      (int*) NULL;
   sol->nweights      = (COLORNWT*) NULL;
   sol->sort_work     =      (int*) NULL;
   sol->grdy_nweights = (COLORNWT*) NULL;
   sol->work_marker   =      (int*) NULL;
   sol->work_path     =      (int*) NULL;
   sol->nodestack     =      (int*) NULL;
   sol->heapref       =      (int*) NULL;
   sol->heap          = (COLORNWTHeap*) NULL;

   COLORclasses_init(&(sol->cclasses));

   sol->lp            = (COLORlp *) NULL;
   sol->dbl_nweights  = (double*) NULL;

   sol->ncount    = 0;
   sol->solcount  = 0;
   sol->freecount = 0;
}

static void print_soldata(soldata* sol)
{
   if (COLORdbg_lvl() > 1) {

      int i;
      printf("WEIGHTS   ");
      for (i = 0;i < sol->ncount;++i) {
         printf(" %5.3g",
                (double) sol->nweights[i] / (double) sol->cutoff);
      }
      printf("\n");

      printf("SOLDATA  ");
      for (i = 0;i < sol->solcount;++i) {
         printf(" %5d",sol->nperm[i]);
      }
      printf("\n");

      printf("TIGHTNESS ");
      for (i = 0;i < sol->ncount;++i) {
         printf(" %5d",sol->tightness[i]);
      }
      printf("\n");
   }
}

static int add_zero_weigthed(soldata* sol)
{
   int i;
   int changes = 0;
   int nnodes = sol->freecount;

   memcpy(sol->sort_work,sol->nperm + sol->solcount,
          sol->freecount * sizeof(int));
   if (COLORdbg_lvl() > 1) {
      printf("Number of free vertices: %d\n",nnodes);
   }
   for (i = 0; i < nnodes; ++i) {
      int v = sol->sort_work[i];
      if (sol->nweights[v] != 0) {
         printf ("Failed to add vertex %d with weight %f upfront.\n",
                 v,COLORunsafe_dbl(sol->nweights[v],sol->cutoff));
      }
      changes += add_iff_free(sol,v);
   }
   return changes;
}

static int remove_zero_weigthed(soldata* sol)
{
   int i;
   int changes = 0;
   int nnodes = sol->solcount;

   memcpy(sol->sort_work,sol->nperm,
          sol->solcount * sizeof(int));
   if (COLORdbg_lvl() > 1) {
      printf("There are %d soldata vertices.\n",nnodes);
   }
   for (i = 0; i < nnodes; ++i) {
      int v = sol->sort_work[i];
      if (sol->nweights[v] == 0.0) {
         if (COLORdbg_lvl() > 1) {
            printf ("Removing zero weight vertex %d from soldata.\n",v);
         }

         changes++;
         remove_from_soldata(sol,v);
      }
   }
   return changes;
}


static int build_LP(soldata*  sol)
{
   int rval = 0;
   int i;
   int ncount = sol->ncount;
   node* nodelist = sol->G->nodelist;

   sol->dbl_nweights = (double*) COLOR_SAFE_MALLOC(ncount,double);
   COLORcheck_NULL(sol->dbl_nweights,"Failed to allocate sol->dbl_nweights");

   COLORlp_free (&(sol->lp));
   rval = COLORlp_init (&(sol->lp), "GreedyLP");
   COLORcheck_rval(rval,"Failed in COLORlp_init for \"GreedyLP\"");

   rval =  COLORlp_objective_sense (sol->lp, -1);
   COLORcheck_rval(rval,"Failed in COLORlp_objective_sense");

   for (i = 0; i < ncount; i++) {
      double w = 0.0;
      rval = COLORlp_addcol (sol->lp, 0, (int *) NULL, (double *) NULL,
                             w, 0.0, 1.0, COLORlp_CONTINUOUS, NULL);
      COLORcheck_rval (rval, "COLORlp_addcol failed");
   }

   for (i = 0; i < ncount; i++) {
      int numnz     = 2;
      int v[2];
      double coef[2] = {1.0, 1.0};
      char sense = COLORlp_LESS_EQUAL;
      int j;
      v[0] = i;
      for( j = 0; j < nodelist[i].degree; ++j) {
         int w_j = nodelist[i].adj[j];
         if (w_j > i) {
            j = nodelist[i].degree;
         } else {
            v[1] = w_j;
            rval = COLORlp_addrow (sol->lp, numnz,v,coef,sense,1,(char *) NULL);
            COLORcheck_rval (rval, "COLORlp_addrow failed");
         }
      }
   }
 CLEANUP:
   return rval;
}


MAYBE_UNUSED static void init_greedy_neighbour_weigths(soldata* sol)
{
   int i;

   node* nodelist = sol->G->nodelist;

   for (i = sol->solcount; i < sol->solcount + sol->freecount;++i) {
      int x = sol->nperm[i];
      int j;
      sol->grdy_nweights[x] = - sol->nweights[x];
      for(j = 0; j < nodelist[x].degree; ++j) {
         int w_j = nodelist[x].adj[j];
         if (is_free(sol,w_j)) {
            sol->grdy_nweights[x] +=
               sol->nweights[w_j];
         }
      }
   }
}

MAYBE_UNUSED static int greedy_improvement_1(soldata* sol, int s)
{
   int i;
   int changes = 0;
   int nnodes = sol->freecount;

   if (!nnodes) {return changes;}

   memcpy(sol->sort_work,sol->nperm + sol->solcount,
          sol->freecount * sizeof(int));

   perm_nwt_rquicksort(sol->sort_work,
                       sol->nweights,
                       sol->freecount);

   for (i = 0; i < nnodes + s; ++i) {
      int v = sol->sort_work[(i + s) % nnodes];
      if (sol->nweights[v] > 0) {
         if(add_iff_free(sol,v)) {
            changes ++;
         }
      }
   }
   return changes;
}

MAYBE_UNUSED static int greedy_improvement_2(soldata* sol, int s)
{
   int i;
   int changes = 0;
   int nnodes = sol->freecount;

   init_greedy_neighbour_weigths(sol);

   memcpy(sol->sort_work,sol->nperm + sol->solcount,
          sol->freecount * sizeof(int));

   perm_nwt_rquicksort(sol->sort_work,
                       sol->grdy_nweights,
                       sol->freecount);


   for (i = 0; i < nnodes; ++i) {
      int v = sol->sort_work[(i + s) % nnodes];
      if (sol->nweights[v] > 0) {
         changes += add_iff_free(sol,v);
      }
   }
   return changes;
}

MAYBE_UNUSED static void adjust_neighbors(soldata* sol, int x)
{
   int i;
   node* nodelist = sol->G->nodelist;
   COLORNWT weight = sol->nweights[x];
   if (weight == 0) {return;}
   for(i = 0; i < nodelist[x].degree; ++i) {
      int w_i = nodelist[x].adj[i];
      if (is_free(sol,w_i)) {
         sol->grdy_nweights[w_i] -= weight;
         COLORNWTheap_decrease_key(sol->heap,
                                   sol->heapref[w_i],
                                   sol->grdy_nweights[w_i]);
      }
   }
}

MAYBE_UNUSED static int greedy_improvement_dyn(soldata* sol, int s)
{
   int  i;
   int  changes = 0;
   int* minv_ptr;

   init_greedy_neighbour_weigths(sol);
   COLORNWTheap_reset(sol->heap);
   for (i = 0; i < sol->freecount;++i) {
      int      v   = sol->nperm[i];
      COLORNWT key = sol->grdy_nweights[v];
      int*     pos = &(sol->heapref[v]);
      minv_ptr     = &(sol->inperm[v]);
/*       printf("Inserting key %lld ptr %p node %d\n", */
/*              (long long)key, (void*) minv_ptr, v); */

      COLORNWTheap_insert(sol->heap, pos,
                          key, (void*) minv_ptr);
   }

   minv_ptr = (int*) COLORNWTheap_min(sol->heap);
   while(minv_ptr) {
      if (!s) {
         int v = sol->nperm[*minv_ptr];
/*          printf("Retrieved key %lld ptr %p node %d\n", */
/*                 (long long)sol->grdy_nweights[v], (void*)minv_ptr, v); */
         if (sol->nweights[v] > 0) {
            int added = add_iff_free(sol,v);
            if (added) {
               changes += added;
               adjust_neighbors(sol, v);
            }
         }
      } else {
         --s;
      }
      minv_ptr = (int*) COLORNWTheap_min(sol->heap);
   }
/*    printf("Retrieved last key\n"); */

   changes += greedy_improvement_2(sol,0);
   return changes;
}


MAYBE_UNUSED static int greedy_improvement_LP(soldata*  sol, int s)
{
   int rval = 0;
   int i;
   COLORNWT scalef;
   int ncount = sol->ncount;
   int nnodes = sol->freecount;
   int* minv_ptr;

   int changes = 0;

   memcpy(sol->sort_work,sol->nperm + sol->solcount,
          sol->freecount * sizeof(int));

   rval = build_LP(sol);
   COLORcheck_rval(rval,"Failed in build_LP");

   COLOR_COLORNWT2double(sol->dbl_nweights,sol->nweights,
                         sol->cutoff, ncount);
   rval = COLORlp_change_objective(sol->lp,0,sol->ncount,sol->dbl_nweights);
   COLORcheck_rval(rval,"Failed in COLORlp_change_objective");

   COLORlp_write(sol->lp,"GreedyLP.lp");

   rval = COLORlp_optimize (sol->lp);
   COLORcheck_rval(rval,"Failed in COLORlp_optimize");

   rval = COLORlp_x (sol->lp, sol->dbl_nweights);
   COLOR_double2COLORNWT(sol->grdy_nweights,&scalef,
                         sol->dbl_nweights,ncount);

   COLORNWTheap_reset(sol->heap);
   for (i = 0; i < sol->freecount;++i) {
      int      v   = sol->nperm[i];
      COLORNWT key = sol->grdy_nweights[v];
      int*     pos = &(sol->heapref[v]);
      minv_ptr     = &(sol->inperm[v]);
/*       printf("Inserting key %lld ptr %p node %d\n", */
/*              (long long)key, (void*) minv_ptr, v); */

      COLORNWTheap_insert(sol->heap, pos,
                          key, (void*) minv_ptr);
   }


   minv_ptr = (int*) COLORNWTheap_min(sol->heap);
   while(minv_ptr) {
      if (!s) {
         int v = sol->nperm[*minv_ptr];
         if (sol->nweights[v] > 0) {
            int added = add_iff_free(sol,v);
            if (added) {
               int cind[1];
               double coef[1] = {1.0};
               int k;
               char sense = COLORlp_EQUAL;
               cind[0] = 1;
               changes++;
               rval = COLORlp_addrow (sol->lp, 1,cind,coef,sense,1,(char *) NULL);
               COLORcheck_rval (rval, "COLORlp_addrow failed");

               rval = COLORlp_optimize (sol->lp);
               COLORcheck_rval(rval,"Failed in COLORlp_optimize");

               rval = COLORlp_x (sol->lp, sol->dbl_nweights);
               COLOR_double2COLORNWT(sol->grdy_nweights,&scalef,
                                     sol->dbl_nweights,ncount);
               for(k = 0; k < nnodes; ++k) {
                  int w = sol->sort_work[k];
                  if (is_free(sol,w)) {
                     COLORNWTheap_relabel(sol->heap,
                                          sol->heapref[w],
                                          sol->grdy_nweights[w]);
                  }
               }
            }
         }
      } else {
         --s;
      }
      minv_ptr = (int*) COLORNWTheap_min(sol->heap);
   }

 CLEANUP:
   return rval;
}


static void reinit_soldatas(soldata*  sol)
{
   int i;
   int ncount = sol->ncount;
   sol->solcount  = 0;
   sol->freecount = ncount;

   for (i = 0;i < ncount;++i) {
      sol->nperm[i] = i;
      sol->inperm[i] = i;
      sol->tightness[i] = 0;
   }
}

static int init_mwis_grdy(soldata*  sol,
                          graph*     G,
                          int          ncount,
                          const COLORNWT nweights[],
                          COLORNWT       cutoff)
{
   int rval = 0;

   sol->ncount    = ncount;
   sol->G         = G;
   sol->nweights  = nweights;
   sol->cutoff    = cutoff;


   if (!sol->nperm) {
      sol->nperm = (int*) COLOR_SAFE_MALLOC(ncount,int);
      COLORcheck_NULL(sol->nperm,"Allocating sol->nperm failed");
   }

   if (!sol->inperm) {
      sol->inperm = (int*) COLOR_SAFE_MALLOC(ncount,int);
      COLORcheck_NULL(sol->inperm,"Allocating sol->inperm failed");
   }

   if (!sol->tightness) {
      sol->tightness = (int*) COLOR_SAFE_MALLOC(ncount,int);
      COLORcheck_NULL(sol->tightness,"Allocating sol->tightness failed");
   }

   if(!sol->sort_work) {
      sol->sort_work = (int*) COLOR_SAFE_MALLOC(ncount,int);
      COLORcheck_NULL(sol->sort_work,"Allocating sol->sort_work failed");
   }

   if (!sol->grdy_nweights) {
      sol->grdy_nweights = (COLORNWT*) COLOR_SAFE_MALLOC(ncount,COLORNWT);
      COLORcheck_NULL(sol->grdy_nweights,"Allocating sol->grdy_nweights failed");
   }
   if (!sol->work_marker) {
      sol->work_marker = (int*) COLOR_SAFE_MALLOC(ncount,int);
      COLORcheck_NULL(sol->work_marker,"Allocating sol->work_marker failed");
   }
   if(!sol->work_path) {
      sol->work_path = (int*) COLOR_SAFE_MALLOC(ncount,int);
      COLORcheck_NULL(sol->work_path,"Allocating sol->work_path failed");
   }
   if(!sol->nodestack) {
      sol->nodestack = (int*) COLOR_SAFE_MALLOC(ncount,int);
      COLORcheck_NULL(sol->nodestack,"Allocating sol->nodestack failed");
   }

   if(!sol->heapref) {
      sol->heapref = (int*) COLOR_SAFE_MALLOC(ncount,int);
      COLORcheck_NULL(sol->heapref,"Allocating sol->heapref failed");
   }

   if(!sol->heap) {
      rval = COLORNWTheap_init(&(sol->heap),ncount);
      COLORcheck_rval(rval, "Failed in COLORNWTheap_init");
   }


   sol->cclasses.cnt = 0;
   sol->greedy_improvement = greedy_improvement_2;

   reinit_soldatas(sol);

/*    greedy_improvement(sol); */

/*    if (COLORdbg_lvl()) { */
/*       print_soldata(sol); */
/*       if (greedy_improvement(sol)) { */
/*          if (COLORdbg_lvl()) printf("ERROR greedy could improve incrementally!\n"); */
/*          print_soldata(sol); */
/*       } */
/*    } */



 CLEANUP:
   if (rval) {
      clean_soldata(sol);
   }
   return rval;
}

static int is_onetight(soldata* sol, int v)
{
   return (sol->tightness[v] == 1);
}

static int is_free_or_onetight(soldata* sol, int v)
{
   return is_onetight(sol,v) || is_free(sol,v);
}

static void remove_from_soldata(soldata* sol, int v)
{
   assert(is_basis(sol, v));

   relax_neighbours(sol,v);
   swap_nodelist(sol,v,sol->nperm[--sol->solcount]);
   sol->freecount++;
   /** v is free now. */
}

static void add_to_soldata(soldata* sol, int v)
{
   assert(is_free(sol,v));

   tighten_neighbours(sol,v);
   swap_nodelist(sol,v,sol->nperm[sol->solcount++]);
   sol->freecount--;
}


static int perform_2_improvement(soldata* sol, int x)
{
   int i;
   int ntight = 0;
   int best_v = -1;
   int best_w = -1;
   const graph* G = sol->G;
   node* nodelist = G->nodelist;
   const COLORNWT* nweights = sol->nweights;
   double best_weight = nweights[x];

   for (i = 0; i < nodelist[x].degree;++i) {
      int v = nodelist[x].adj[i];

      if(is_free_or_onetight(sol,v)) {
         ntight++;
         if (ntight > 1) {
            /* Now find tight non-neighbor w:*/
            int l_i,n_i = 0;
            /* Loop through one-tight neighbors of x.*/
            for(l_i = 0; 
                (l_i < nodelist[x].degree )  && (n_i < nodelist[v].degree);
                l_i++) 
            {
               int w_l = nodelist[x].adj[l_i];
               if (w_l == v) {continue;}

               if (is_free_or_onetight(sol,w_l)) {
                  int w_n = nodelist[v].adj[n_i];

                  /* Look wether w_n occurs in neighborhood of v.*/
                  while (w_n < w_l && n_i < nodelist[v].degree) {
                     n_i++;
                     if (n_i < nodelist[v].degree) {
                        w_n = nodelist[v].adj[n_i];
                     }
                  }
                  if (w_n > w_l || n_i == nodelist[v].degree ) {
                     /* feasible w found.*/
                     double swap_weight = nweights[v] + nweights[w_l];
                     if (swap_weight > best_weight) {
                        best_v = v;
                        best_w = w_l;
                        best_weight = swap_weight;
                     }
                  }
               }
            }
         }
      }
   }
   if (best_v != -1) {

      if(COLORdbg_lvl() > 1) printf("Found improving swap %d <- (%d,%d) (delta %f)\n",x,
                                    best_v,best_w,best_weight - nweights[x]);

      swap_nodelist(sol,x, best_v);

      relax_neighbours(sol,x);
      tighten_neighbours(sol,best_v);

      tighten_neighbours(sol,best_w);
      swap_nodelist(sol,best_w,sol->nperm[sol->solcount]);
      sol->solcount++;
      sol->freecount--;

      if(COLORdbg_lvl() > 1) print_soldata(sol);

      return 1;
   }
   return 0;
}


static int perform_2_improvements(soldata* sol)
{
   int i;
   int totnswaps = 0;
   int iternswaps;
   if (COLORdbg_lvl() > 1) {
      printf("Starting perform_2_improvements.\n");
   }
   do {
      iternswaps = 0;

      for (i = 0; i < sol->solcount;++i) {
         int changed = perform_2_improvement(sol,sol->nperm[i]);
         if (changed) {
            --i;
            iternswaps += changed;
         }
      }
      if (iternswaps) {
         int changed = sol->greedy_improvement(sol,0);
         totnswaps += iternswaps;
         if (changed) {
            if(COLORdbg_lvl() > 1) print_soldata(sol);
         }
      }
   }  while (iternswaps);
   return totnswaps;
}

static void mark_neighbors(int work_marker[], const node* n)
{
   int e_i;
   for (e_i = 0; e_i < n->degree; ++e_i) {
      int w = n->adj[e_i];
      (work_marker[w])++;
   }
}


static int perform_1_2_path(soldata* sol, int* nodestack, int v)
{
   int i;
   int rm_i = 0;
   int add_i = sol->ncount;
   int s_i = 0;
   COLORNWT weight  = 0;
   const graph* G     = sol->G;
   const node*  nodes = G->nodelist;
   int changes = 0;


   if (COLORdbg_lvl() > 1) {
      printf("Starting 1_2 path search with %d\n",v);
   }

   weight = 0;
   for (i = 0; i < sol->ncount; ++i) {
      sol->work_marker[i] = 0;
      sol->work_path[i] = -1;
      nodestack[i] = -1;
   }

   sol->work_path[--add_i] = v;
   weight += sol->nweights[v];
   mark_neighbors(sol->work_marker,&(nodes[v]));
   nodestack[s_i++] = v;

   while (s_i) {
      int x = nodestack[--s_i];
      nodestack[s_i] = -1;
      if(is_basis(sol,x)) {
         int e_i;
         for( e_i = 0; e_i < nodes[x].degree; ++e_i) {
            int y = nodes[x].adj[e_i];
            if (is_onetight(sol,y) && sol->nweights[y] > 0) {
               if (! sol->work_marker[y]) {
                  sol->work_path[--add_i] = y;
                  weight += sol->nweights[y];
                  mark_neighbors(sol->work_marker,&(nodes[y]));
               }
            }
         }
      } else {
         int e_i;
         for( e_i = 0; e_i < nodes[x].degree; ++e_i) {
            int y = nodes[x].adj[e_i];
            if (is_basis(sol,y)) {
               sol->work_path[rm_i++] = y;
               nodestack[s_i++] = y;
               weight -= sol->nweights[y];
            }
         }
      }
   }
   if (weight > 0) {
      if (COLORdbg_lvl() > 1) {
         printf("Found replacement (");
         for (rm_i = 0; sol->work_path[rm_i] != -1;++rm_i) {
            printf(" %d", sol->work_path[rm_i]);
         }
         printf(") <- (");
         for (add_i = sol->ncount-1; sol->work_path[add_i] != -1;--add_i) {
            printf(" %d",sol->work_path[add_i]);
         }
         printf(") weight %lld\n",(long long) weight);
      }

      for (rm_i = 0; sol->work_path[rm_i] != -1;++rm_i) {
         int x = sol->work_path[rm_i];
         remove_from_soldata(sol,x);
      }
      for (add_i = sol->ncount-1; sol->work_path[add_i] != -1;--add_i) {
         int x = sol->work_path[add_i];
         add_to_soldata(sol,x);
      }
      changes++;
   }
   return changes;
}

static int perform_1_2_paths(soldata* sol)
{
   int rval = 0;
   int i;
   int   nnonfree     = sol->ncount - sol->solcount - sol->freecount;
   int changes    = 0;
   if (!nnonfree) return rval;


   memcpy(sol->sort_work ,
          sol->nperm + sol->solcount + sol->freecount,
          nnonfree*sizeof(int));


   for (i = 0; i < nnonfree;++i) {
      if (sol->tightness[sol->sort_work[i]] == 2) {
         sol->grdy_nweights[sol->sort_work[i]] = sol->nweights[sol->sort_work[i]];
      } else {
         sol->grdy_nweights[sol->sort_work[i]] = .0;
      }
   }
   perm_nwt_rquicksort (sol->sort_work, sol->grdy_nweights, nnonfree);

   for (i = 0; i < nnonfree;++i) {
      int v = sol->sort_work[i];

      if (sol->tightness[v] == 2 && sol->nweights[v] > DBL_EPSILON) {
         changes += perform_1_2_path(sol,sol->nodestack,v);
      }
   }
   return changes;
}

static COLORNWT soldata_value(soldata* sol)
{
   int i;
   COLORNWT sval = 0;
   for (i = 0; i < sol->solcount;++i)
      {
         sval += sol->nweights[sol->nperm[i]];
   }

   if (COLORdbg_lvl() > 1)
   {
      printf("mwis_greedy found  %20.15f ( %lld / %lld ).\n",
             COLORunsafe_dbl(sval,sol->cutoff),
             (long long) sval, (long long)sol->cutoff);
   }

   return sval;
}

int COLORcheck_set(COLORset* set, int ncount, int ecount, const int elist[])
{
   int rval = 0;
   int i;
   int* coloring = (int*) COLOR_SAFE_MALLOC(ncount,int);
   COLORcheck_NULL(coloring,"Could not allocate *newsets");

   for (i = 0; i < ncount;++i) {
      coloring[i] = 0;
   }

   for (i = 0; i < set->count;++i) {
      coloring[set->members[i]] = 1;
   }

   for (i = 0; i < ecount; ++i) {
      if (coloring[elist[2*i]] == 1 && coloring[elist[2*i+1]] ==1) {
         fflush(stdout);
         fprintf(stderr,"ILLEGAL COLORING FOUND!\n");
         rval++;
      }
   }
 CLEANUP:
   if (coloring) free(coloring);
   return rval;
}

static int vertex_comparator(const void* v1,const void* v2)
{
   int i1 = *(const int*) v1;
   int i2 = *(const int*) v2;

   return i1-i2;
}

static int add_soldata(soldata* sol)
{
   int rval = 0;
   int* newmembers = (int*) NULL;
   COLORclasses* cclasses = &(sol->cclasses);

   if (cclasses->cnt+1 > cclasses->allocnt) {
      rval = COLORclasses_expand(cclasses);
      COLORcheck_rval(rval,"Failed in COLORclasses_expand");
   }

   cclasses->sets[cclasses->cnt].count = sol->solcount;
   newmembers  =
      (int*) COLOR_SAFE_MALLOC(sol->solcount,int);
   COLORcheck_NULL(newmembers,
                   "Failed to allocate classes->sets[classes->cnt].members");
   memcpy(newmembers,
          sol->nperm,sol->solcount*sizeof(int));

   qsort(newmembers,sol->solcount,sizeof(int),vertex_comparator);

   cclasses->sets[cclasses->cnt].members = newmembers;
   ++(cclasses->cnt);
   {
      int i;
      printf("NEW SET ");
      for (i = 0; i < sol->solcount;++i) {
         printf(" %d",newmembers[i]);
      }
      printf("\n");
   }
 CLEANUP:
   if (rval) {
      if (newmembers) {free(newmembers);}
   }
   return rval;
}


static int transfer_soldata(COLORset** newsets,
                             int*       nnewsets,
                             soldata*  sol)
{
   int rval = 0;
   COLORclasses* cclasses = &(sol->cclasses);

   if (cclasses->cnt == 0) goto CLEANUP;

   if (COLORdbg_lvl() > 1) {
      printf("Transferring %d greedy soldatas.\n",cclasses->cnt);
   }

   *nnewsets = cclasses->cnt;
   *newsets = (COLORset *) COLOR_SAFE_MALLOC(*nnewsets,COLORset);
   COLORcheck_NULL(*newsets,"Could not allocate *newsets");

   memcpy(*newsets,cclasses->sets,*nnewsets * sizeof(COLORset));
   COLORcheck_NULL(*newsets,"Could not allocate *newsets");

   COLORclasses_reset_without_free(cclasses);

CLEANUP:
   if (rval) {
      if (*newsets) {
         if ((*newsets)[0].members) { free((*newsets)[0].members);}
         free(*newsets); *newsets = (COLORset*) NULL;
      }
   }
   return rval;
}

static int inspect_soldata(soldata* sol, COLORNWT* best_sval,
                            const char* logstring)
{
   int    rval = 0;
   COLORNWT sval = soldata_value(sol);
   if (COLORdbg_lvl() > 1) printf("%s: %lld\n",logstring,(long long) sval);

   if(sval > sol->cutoff && sval > *best_sval) {
      add_zero_weigthed(sol);
      rval = add_soldata(sol);
      remove_zero_weigthed(sol);

      COLORcheck_rval(rval,"Failed in add_soldata");
      *best_sval = sval;
   } else if (sval > *best_sval) {
      *best_sval = sval;
      print_soldata(sol);
   }
 CLEANUP:
   return rval;
}



static int repeated_greedy_followed_by_ls(soldata*  sol)
{
   int rval = 0;
   COLORNWT best_sval = 0;
   int changes;
   int start_vertex;
   int nsoldatas = 0;
   int last_improving_start = 0;
   int last_valid_start = -1;
   int num_starts = sol->ncount; /* (int) sqrt((double) sol->ncount); */
   for (start_vertex = 0; start_vertex  < num_starts; ++start_vertex) {
      reinit_soldatas(sol);
      sol->solcount = 0;
      sol->greedy_improvement(sol,start_vertex);
      if ( !(changes = sol->solcount)) continue;

      inspect_soldata(sol,&best_sval,"GREEDY MWIS");

      while (changes) {

         changes = perform_2_improvements(sol);
         COLORcheck_rval(rval,"perform_2_improvements");
         if (changes) {
            inspect_soldata(sol,&best_sval,"GREEDY 2-SWAPS");
         }

         if (best_sval < sol->cutoff) {
            int change = perform_1_2_paths(sol);
            changes += change;
            COLORcheck_rval(rval,"perform_1_2_paths");

            if (change) {
               sol->greedy_improvement(sol,0);
               inspect_soldata(sol,&best_sval,"GREEDY 1-2-SWAPS");
            }
         }
      }

      if (1) {
         int zero_added = add_zero_weigthed(sol);
         if (COLORdbg_lvl() > 1 && zero_added) {
            printf("Added %d zero weighted vertices.\n",zero_added);
         }
      }
      if (sol->cclasses.cnt > nsoldatas) {
         nsoldatas = sol->cclasses.cnt;
         last_improving_start = start_vertex;
         if (last_valid_start == -1)
            last_valid_start = start_vertex;
/*          num_starts /=2; */
         num_starts = sol->ncount / nsoldatas;
      }
   }
   printf("Best greedy:   %13.10e ( %lld / %lld ) , number of greedy soldatas: %d, first valid it. %d last improving iteration %d\n",
          COLORsafe_lower_dbl(best_sval,sol->cutoff),(long long ) best_sval,(long long ) sol->cutoff,
          sol->cclasses.cnt, last_valid_start,last_improving_start);
   fflush(stdout);
 CLEANUP:
   return rval;
}

int COLORstable_LS(MWISls_env** env,
                   COLORset** newsets, int* nnewsets, int ncount,
                   int ecount, const int elist[],
                   const COLORNWT nweights[], COLORNWT cutoff)
{
   int rval = 0;

   graph*     G  = (graph*) NULL;
   soldata* sol = (soldata*) NULL;

   if (! *env) {
      (*env) = (MWISls_env*) COLOR_SAFE_MALLOC (1,MWISls_env);
      COLORcheck_NULL(*env,"Allocating *env failed.");

      set_ptrs_to_zero(&((*env)->sol));

      G   = &((*env)->G);
      G->nodelist = (node*) NULL;
      G->adjspace = (int*) NULL;
      COLORadjgraph_free(G);

      rval =  COLORadjgraph_build(G, ncount, ecount, elist);
      COLORcheck_rval(rval,"COLORbuild_adjgraph failed");

      rval = COLORadjgraph_simplify(G);

      COLORadjgraph_sort_adjlists_by_id(G);
   }

   G   = &((*env)->G);
   sol = &((*env)->sol);


   rval = init_mwis_grdy(sol,G,ncount,nweights,cutoff);
   COLORcheck_rval(rval,"init_mwis_grdy");

   if (COLORdbg_lvl() > 1)
      printf("Starting new round of repeated_greedy_followed_by_ls...\n");

   sol->greedy_improvement = greedy_improvement_2;
   repeated_greedy_followed_by_ls(sol);

   if(!sol->cclasses.cnt) {
      sol->greedy_improvement = greedy_improvement_1;
      repeated_greedy_followed_by_ls(sol);
   }

   if(!sol->cclasses.cnt) {
      sol->greedy_improvement = greedy_improvement_dyn;
      repeated_greedy_followed_by_ls(sol);
   }

   /** In fact the LP-based greedy & local is essentially what
       gurobi does (reduced to a predetermined set of branches.
       But without cutting and bounding  it's not faster than gurobi.
       Instead we might try  a fixed number of DFS branches within
       stable.c.
   */
/*    if (!sol->cclasses.cnt) { */
/*       sol->greedy_improvement = greedy_improvement_LP; */
/*       repeated_greedy_followed_by_ls(sol); */
/*    } */


   if(sol->cclasses.cnt) {
      transfer_soldata(newsets,nnewsets,sol);
      rval = COLORcheck_set(newsets[0],ncount,ecount,elist);
      COLORcheck_rval(rval,"COLORcheck_set failed");
   }


 CLEANUP:
   if (rval) {
      COLORfree_sets (newsets,nnewsets);
      clean_soldata(sol);
      COLORadjgraph_free(G);
   }

   return rval;
}

static int num_weighted_neighbors(soldata* sol, COLORNWT* nweights, int v)
{
   int j;
   int numwn = 0;
   node* nodelist = sol->G->nodelist;

   for (j = 0;j < nodelist[v].degree; ++j) {
      int w_j = nodelist[v].adj[j];
      if (nweights[w_j] > 0) {
         numwn++;
      }
   }
   return numwn;
}

static int decrease_nwn_of_neighbors(soldata* sol, int v)
{
   int rval;
   int k;
   
   node* nodelist = sol->G->nodelist;

   for (k = 0;k < nodelist[v].degree; ++k) {
      int w_k = nodelist[v].adj[k];
      if (sol->heapref[w_k] != -1) {
/*          printf("Decreasing key of node %d to %d.\n", */
/*          w_k, sol->work_marker[w_k] - 1); */
         --(sol->work_marker[w_k]);
         rval = COLORNWTheap_decrease_key(sol->heap,
                                          sol->heapref[w_k],
                                          sol->work_marker[w_k]);
         COLORcheck_rval(rval,"Failed in COLORNWTheap_decrease_key");
      }
   }
 CLEANUP:
   return rval;
}

int COLORstable_round_down_weights(MWISls_env* env,
                                   COLORNWT nweights[],
                                   COLORNWT cutoff)
{
   int rval = 0;
   int i;
   int ncount = env->G.ncount;
   const node* nodelist    = env->G.nodelist;
   soldata* sol = &(env->sol);
   COLORNWT lb = 0;
   COLORNWT lb_frac = 0;
/*    COLORNWT min_frac = (COLORNWT) ncount; */
   COLORNWT min_frac = (COLORNWT) 1; 

   int* minv_ptr;

   for (i = 0; i < ncount; ++i) {
      lb += nweights[i];
      sol->heapref[i] = -1;
      sol->work_marker[i] = 0; // using work_marker to store keys here;
   }

   /** The remainder of lb/cutoff (decremented by one) can be used to
       round down the node weihts without weakening the lower bound.
       E.g. a fractional lower bound of 15.43 will be decreased to 15
       + \epsilon while maintaining the integral lower bound of 16.
       This way we may spare some column generation iterations.
   */

   lb_frac = lb % cutoff  - min_frac;

   COLORNWTheap_reset(sol->heap);
   for (i = 0; (i < ncount) && (lb_frac > 0); ++i) {
      int*     pos = &(sol->heapref[i]);
      if (nweights[i] > 0) {
         sol->work_marker[i] = 1;
         sol->work_marker[i] += num_weighted_neighbors(sol,nweights,i);
         minv_ptr     = &(sol->inperm[i]);
         rval  = COLORNWTheap_insert(sol->heap, pos,
                                     sol->work_marker[i], 
                                     (void*) minv_ptr);
/*          printf("Inserted node %d with key %d at pos %d, %d.\n", */
/*                 i,sol->work_marker[i],*pos,sol->heapref[i]); */
         COLORcheck_rval(rval,"Failed in COLORNWTheap_insert");
      }
   }

   minv_ptr = (int*) COLORNWTheap_min(sol->heap);
   while (lb_frac > 0 && minv_ptr) {
      int j;
      int v = sol->nperm[*minv_ptr];
      COLORNWT adjust;
      int start = 0, end = 0;
      /* At any time there must not be nodes v with nweights[v] == 0
         in the heap!
      */
      assert (nweights[v] > 0);

      sol->heapref[v] = -1;

      sol->sort_work[end++] = v;
      for (j = 0; j < nodelist[v].degree; ++j) {
         int w_j = nodelist[v].adj[j];
         if (nweights[w_j] > 0) {
            sol->sort_work[end++] = w_j;
         }
      }
      perm_nwt_rquicksort(sol->sort_work,nweights,end);

      while ( (lb_frac > 0) && 
              (start < end) &&
              ((adjust = nweights[sol->sort_work[end-1]]) > 0) ) {
         adjust = COLORNWTmin(adjust,lb_frac/(end-start));

         lb_frac -= adjust * (end-start);

         for (j = start; j < end; ++j) {
            int w_j= sol->sort_work[j];
            assert (nweights[w_j] >= adjust);
            nweights[w_j] -= adjust;
            if (!nweights[w_j]) {
               if (sol->heapref[w_j] != -1) {
                  /*                   printf("Removing node %d from pos %d.\n",w_j,sol->heapref[w_j]); */
                  rval = COLORNWTheap_remove(sol->heap,sol->heapref[w_j]);
                  COLORcheck_rval(rval,"Failed in COLORNWTheap_remove");
                  sol->heapref[w_j] = -1;
               }
               decrease_nwn_of_neighbors(sol,w_j);
            }
         }
         if (sol->sort_work[end-1] == v ) {
            start = end;
         } else {
            end--;
         }
         if ( lb_frac <= (end-start)) {
            start = end - lb_frac;
         }
      }
      minv_ptr = (int*) COLORNWTheap_min(sol->heap);
   }
 CLEANUP:
   return rval;
}

int COLORstable_free_ls_env(MWISls_env** env)
{
   if (*env) {
      COLORadjgraph_free(&((*env)->G));
      clean_soldata(&((*env)->sol));
      free(*env);
      *env = (MWISls_env*) NULL;
   }
   return 0;
}
