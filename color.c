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
#include <sys/resource.h>
#include <time.h>
#include <assert.h>
#include <limits.h>

#include "heap.h"
#include "lp.h"
#include "graph.h"
#include "color.h"

#include "mwis.h"
#include "plotting.h"

static char *edgefile = (char *) NULL;
static char *outfile = (char *) NULL;
static char *cclasses_outfile = (char *) NULL;
static char *cclasses_infile = (char *) NULL;
static char *color_infile = (char *) NULL;


static int debug = 0;
static int write_mwis = 0;
static int ncolordata = 0;
static int initial_upper_bound = INT_MAX;

typedef struct colordata colordata;
struct colordata {

   /* The id of the node in the B&B tree.*/
   int id;
   int depth;

   enum {
      initialized       = 0,
      LP_bound_computed = 1,
      finished          = 2,
   } status;

   /* The instance graph */
   int ncount;
   int ecount;
   int *elist;
   int *orig_node_ids;


   /* The column generation LP. */
   COLORlp * lp;
   double *coef;
   double *pi;

   COLORNWT  mwis_pi_scalef;
   COLORNWT *mwis_pi;
   COLORNWT  lower_bound;
   COLORNWT  upper_bound;
   COLORNWT  lower_scaled_bound;
   double    dbl_lower_bound;

   /* The MWIS instances. */
   MWISenv*  mwis_env;
   COLORset *cclasses;
   int       ccount;
   int       gallocated;
   COLORset *newsets;
   int       nnewsets;

   COLORset *bestcolors;
   int       nbestcolors;


   const COLORset *debugcolors;
   int             ndebugcolors;
   int   opt_track;

   int maxiterations;
   int retirementage;

   int v1,v2;
   colordata*       parent;
   colordata*       same_child;
   colordata*       diff_child;

   char  pname[256];
};



int main (int ac, char **av);
static int parseargs (int ac, char **av);
static void usage (char *f);
static int get_problem_name(char* pname,const char* efname);
static int build_lp(colordata* cd);
static int compute_lower_bound(colordata* cd);
static int compute_coloring(colordata* cd);
static int grab_integral_solution(colordata* cd, double* x, double tolerance);

static void init_colordata(colordata* cd)
{
   cd->id = ncolordata++;
   cd->depth = 0;
   cd->status = initialized;
   sprintf(cd->pname,"colorprob_%d",cd->id);
   cd->ncount = 0;
   cd->ecount = 0;
   cd->elist = (int *) NULL;
   cd->orig_node_ids = (int*) NULL;

   cd->coef = (double *) NULL;
   cd->pi = (double *) NULL;
   cd->mwis_pi_scalef = 1.0;
   cd->mwis_pi = (COLORNWT *) NULL;
   cd->upper_bound = initial_upper_bound;
   cd->lower_bound = 1;
   cd->dbl_lower_bound = 1.0;
   cd->lower_scaled_bound = 1;


   cd->lp       = (COLORlp *) NULL;
   cd->mwis_env  = (MWISenv*) NULL;

   cd->newsets = (COLORset*) NULL;
   cd->nnewsets = 0;

   cd->ccount = 0;
   cd->cclasses = (COLORset *) NULL;
   cd->gallocated = 0;

   cd->bestcolors = (COLORset *) NULL;
   cd->nbestcolors = 0;

   cd->debugcolors  = (COLORset *) NULL;
   cd->ndebugcolors = 0;
   cd->opt_track    = 0;

   cd->maxiterations = 1000000;
   cd->retirementage = 1000000;

  cd->v1 = cd->v2 = -1;
   cd->parent = (colordata*) NULL;

   cd->same_child = (colordata*)NULL;
   cd->diff_child = (colordata*)NULL;
}

static void free_lbcolordata(colordata* cd)
{
   if (cd->lp) COLORlp_free (&(cd->lp));
   if (cd->coef) {
      free (cd->coef);
      cd->coef = (double*) NULL;
   }
   if (cd->pi) {
      free (cd->pi);
      cd->pi   = (double*) NULL;
   }
   if (cd->mwis_pi) {
      free(cd->mwis_pi);
      cd->mwis_pi = (COLORNWT*) NULL;
   }
   COLORfree_sets(&(cd->newsets),&(cd->nnewsets));
   COLORfree_sets(&(cd->cclasses),&(cd->gallocated));
   

   COLORstable_freeenv(&(cd->mwis_env));

}

static void free_colordata(colordata* cd)
{
   if (cd->same_child) {
      free_colordata(cd->same_child);
      free(cd->same_child);
   }
   if (cd->diff_child) {
      free_colordata(cd->diff_child);
      free(cd->diff_child);
   }

   free_lbcolordata(cd);
   if (cd->elist) free (cd->elist);
   if (cd->orig_node_ids) free(cd->orig_node_ids);

    COLORfree_sets(&(cd->bestcolors),&(cd->nbestcolors));

}

int COLORdbg_lvl() {
   return debug;
}

static int print_colors(COLORset* cclasses, int ccount)
{
   int i,j;
   for (i = 0; i < ccount; i++) {
      printf ("Color %d:", i);
      for (j = 0; j < cclasses[i].count; j++) {
         printf (" %d",cclasses[i].members[j]);
      }
      printf ("\n");
   }
   fflush (stdout);
   return 0;
}

static void compute_objective(colordata* cd)
{
   int i;

   cd->lower_scaled_bound = .0;
   for (i = 0; i < cd->ncount;++i) {
      cd->lower_scaled_bound += (double) cd->mwis_pi[i];
   }

   cd->dbl_lower_bound = COLORsafe_lower_dbl(cd->lower_scaled_bound,cd->mwis_pi_scalef);
   printf("Current primal LP objective: %f (%lld / %lld).\n",
          cd->dbl_lower_bound,
          (long long) cd->lower_scaled_bound,
          (long long) cd->mwis_pi_scalef );

   cd->lower_bound = ceil(cd->dbl_lower_bound);

}

static void make_pi_feasible(colordata* cd)
{
   int c;
   int current_rounding = fegetround();
   fesetround(FE_UPWARD);

   for (c = 0; c < cd->ccount;++c) {
      int i;
      double colsum = .0;
      double newcolsum = .0;
      for (i = 0; i < cd->cclasses[c].count;++i) {
         if (signbit(cd->pi[cd->cclasses[c].members[i]])) {
            cd->pi[cd->cclasses[c].members[i]] = 0.0;
         }
         colsum += cd->pi[cd->cclasses[c].members[i]];
      }
      if (colsum > 1.0) {
         fesetround(FE_DOWNWARD);

         for (i = 0; i < cd->cclasses[c].count;++i) {
            cd->pi[cd->cclasses[c].members[i]] /= colsum;
            newcolsum += cd->pi[cd->cclasses[c].members[i]];
         }
         if (COLORdbg_lvl()> 1) {printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n",c,colsum,newcolsum);}
         fesetround(FE_UPWARD);
      }
      fesetround(current_rounding);
   }
}


static void reset_ages(COLORset* cclasses, int ccount)
{
   int i;
   for (i = 0; i < ccount; ++i) {
      cclasses[i].age = 0;
   }
}

static int grow_ages(colordata* cd)
{
   int rval = 0;
   int i;
   double* x = (double*) NULL;

   x = (double*) COLOR_SAFE_MALLOC(cd->ccount, double);
   COLORcheck_NULL(x,"Failed to allocate x");

   rval = COLORlp_x(cd->lp,x);
   COLORcheck_rval(rval,"Failed in COLORlp_x");

   for (i = 0; i < cd->ccount; ++i) {
      if (x[i] < DBL_EPSILON) {
         ++(cd->cclasses[i].age);
      } else {
         cd->cclasses[i].age = 0;
      }
   }

 CLEANUP:
   if (x) free(x);
   return rval;
}

static int delete_old_colorclasses(colordata* cd)
{
   int rval   = 0;
   int i;
   int numdel = 0;


   for (i = 0; i < cd->ccount; ++i) {
      if (cd->cclasses[i].age > cd->retirementage) {
         /* swap current class with last and delete it. */
         --(cd->ccount);
         COLORfree_set(&(cd->cclasses[i]));
         memcpy(&(cd->cclasses[i]),&(cd->cclasses[cd->ccount]),sizeof(COLORset));
         COLORinit_set(&(cd->cclasses[cd->ccount]));
         /* Ensure that the formerly last class is considered:*/
         --i;
         numdel++;
      }
   }

   printf("Deleted %d out of %d columns with age > %d. Rebuilding LP from scratch.\n",
          numdel, numdel + cd->ccount, cd->retirementage);

   if (numdel) {
      printf("Deleted %d out of %d columns with age > %d. Rebuilding LP from scratch.\n",
             numdel, numdel + cd->ccount, cd->retirementage);
      COLORlp_free (&(cd->lp));

      rval = build_lp(cd);
      COLORcheck_rval(rval, "Failed in build_lp");
   }


 CLEANUP:

   return rval;
}

static int write_mwis_instances(colordata* cd)
{
   int rval = 0;
   if ( (cd->id == 0) && write_mwis) {
      char   mps_fname[256];

      sprintf(mps_fname,"%s.mwis.lp",cd->pname);
      rval = COLORstable_write_mps(mps_fname,cd->ncount,cd->ecount,cd->elist,
                                   cd->mwis_pi,cd->mwis_pi_scalef);
      COLORcheck_rval(rval,"Failed in COLORstable_write_mps");

      sprintf(mps_fname,"%s.mwis.dimacs",cd->pname);
      rval = COLORstable_write_dimacs(mps_fname,cd->ncount,cd->ecount,cd->elist,
                                      cd->mwis_pi,cd->mwis_pi_scalef);
      COLORcheck_rval(rval,"Failed in COLORstable_write_dimacs");

      sprintf(mps_fname,"%s.mwclq.dimacs",cd->pname);
      rval = COLORstable_write_dimacs_clique(mps_fname,
                                             cd->ncount,cd->ecount,cd->elist,
                                             cd->mwis_pi,cd->mwis_pi_scalef);
      COLORcheck_rval(rval,"Failed in COLORstable_write_dimacs_clique");

/*       sprintf(mps_fname,"%s.graphviz.dot",cd->pname); */
/*       COLORplot_graphviz(mps_fname,cd->ncount,cd->ecount,cd->elist,(int*) NULL); */
/*       COLORcheck_rval(rval,"Failed in COLORplot_graphviz"); */
   }
 CLEANUP:
   return rval;
}

static int write_snapshot(colordata* cd, int add_timestamp)
{
   int rval = 0;
   if (cclasses_outfile != (char*) NULL) {
      char   fname[256];
      const colordata* root_cd = cd;

      while (root_cd->parent) {
         root_cd = root_cd->parent;
      }

      if (add_timestamp) {
         char   timestr[256];
         time_t t;
         t = time(NULL);

         strftime(timestr, 256, "%Y%m%d%H%M", localtime(&t));
         sprintf(fname,"%s.%s",cclasses_outfile,timestr);
      } else {
         sprintf(fname,"%s",cclasses_outfile);
      }
      rval = COLORstable_write_stable_sets(cd->cclasses,cd->ccount,cd->ncount,
                                           fname,root_cd->pname);
      COLORcheck_rval(rval,"Failed in COLORstable_write_stable_sets");
   }
 CLEANUP:
   return rval;
}


static void print_ages(colordata* cd)
{
   int i;
   printf("AGES:");
   for (i = 0; i < cd->ccount; ++i) {
      printf(" %4d",cd->cclasses[i].age);
   }
   printf("\n");

}



static int concat_newsets(colordata* cd)
{
   int rval = 0;
   COLORset* tmpsets = (COLORset*) NULL;
   int i;
   if (cd->nnewsets == 0) return rval;

   reset_ages(cd->newsets,cd->nnewsets) ;

   if (cd->ccount + cd->nnewsets > cd->gallocated) {
      /* Double size */
      tmpsets = COLOR_SAFE_MALLOC(2 * cd->gallocated,COLORset);
      COLORcheck_NULL (tmpsets, "out of memory for tmpsets");
      memcpy(tmpsets,cd->cclasses, cd->ccount * sizeof(COLORset));
      free(cd->cclasses);
      cd->gallocated *= 2;
      cd->cclasses = tmpsets;
      tmpsets = NULL;
   }
   memcpy(cd->cclasses + cd->ccount, cd->newsets, cd->nnewsets * sizeof(COLORset));
   cd->ccount += cd->nnewsets;
   for (i = cd->ccount; i < cd->gallocated; ++i) {
      cd->cclasses[i].count = 0;
      cd->cclasses[i].members = (int*) NULL;
   }
 CLEANUP:
   if (rval) {
      if (cd->cclasses) free(cd->cclasses);
   }
   if (cd->newsets) free(cd->newsets);
   cd->newsets = (COLORset*) NULL;
   cd->nnewsets = 0;
   return rval;
}

static int build_lp(colordata* cd)
{
    int i;

    int rval = COLORlp_init (&(cd->lp), "colorme");
    COLORcheck_rval (rval, "COLORlp_init failed");

    for (i = 0; i < cd->ncount; i++) {
       char sense = 'G';
        rval = COLORlp_addrow (cd->lp, 0, (int *) NULL, (double *) NULL, sense,
                               1.0, (char*) NULL);
        COLORcheck_rval (rval, "COLORlp_addrow failed");
    }

    cd->coef = (double *) realloc (cd->coef,cd->ncount * sizeof (double));
    COLORcheck_NULL (cd->coef, "out of memory for coef");
    for (i = 0; i < cd->ncount; i++) cd->coef[i] = 1.0;

    for (i = 0; i < cd->ccount; i++) {
       rval = COLORlp_addcol (cd->lp, cd->cclasses[i].count, cd->cclasses[i].members,
                              cd->coef, 1.0, 0.0, 1.0, COLORlp_CONTINUOUS, (char*) NULL);
       if (rval) COLORlp_printerrorcode (rval);
       COLORcheck_rval (rval, "COLORlp_addcol failed");
    }

    if (COLORdbg_lvl() > 1) {
       rval = COLORlp_write (cd->lp, "look.lp");
       COLORcheck_rval (rval, "COLORlp_write failed");
    }

    cd->pi = (double *) realloc (cd->pi,cd->ncount * sizeof (double));
    COLORcheck_NULL (cd->pi, "out of memory for pi");
 CLEANUP:
    if (rval) {
       COLORlp_free(&(cd->lp));
       if (cd->coef) free(cd->coef);
       if (cd->pi) free(cd->pi);
    }
    return rval;
}

COLOR_MAYBE_UNUSED static int heur_colors_with_stable_sets(colordata* cd)
{
   int rval = 0;
   double incumbent;
   double* colsol;
   int* colored = (int*) NULL;
   int i;

   colsol = (double*) COLOR_SAFE_MALLOC(cd->ccount,double);
   COLORcheck_NULL(colsol,"Failed to allocate colsol");

   colored = (int*) COLOR_SAFE_MALLOC(cd->ncount,int);
   COLORcheck_NULL(colored,"Failed to allocate colored");

   for (i = 0; i < cd->ncount; ++i) {colored[i] = 0;}

   COLORlp_set_all_coltypes(cd->lp,GRB_BINARY);
   COLORcheck_rval (rval, "COLORlp_set_all_coltypes");

   COLORlp_write (cd->lp, "lpheur.lp");

   rval = COLORlp_setnodelimit(cd->lp,1);
   COLORcheck_rval(rval,"COLORlp_setnodelimit failed");

   rval = COLORlp_optimize(cd->lp);
   COLORcheck_rval (rval, "COLORlp_optimize failed");

   rval =  COLORlp_x (cd->lp, colsol);
   COLORcheck_rval (rval, "COLORlp_x failed");

   rval = COLORlp_objval (cd->lp, &incumbent);
   COLORcheck_rval (rval, "COLORlp_objval failed");

   grab_integral_solution(cd,colsol,COLORlp_int_tolerance());

   COLORlp_set_all_coltypes(cd->lp,GRB_CONTINUOUS);
   COLORcheck_rval (rval, "COLORlp_set_all_coltypes");

   printf ("Found lower bound of %lld and upper bound of %g.\n",
           (long long) cd->lower_bound, incumbent);

   rval = COLORcheck_coloring(cd->bestcolors,cd->nbestcolors,
                              cd->ncount, cd->ecount, cd->elist);
   COLORcheck_rval(rval,"ERROR: An incorrect coloring was created.");

   print_colors(cd->bestcolors,cd->nbestcolors);
 CLEANUP:
   if (colsol) free(colsol);
   if (colored) free(colored);
   return rval;
}

static void COLORset_SWAP(COLORset *c1,COLORset *c2,COLORset *t)
{
   memcpy(t,c2,sizeof(COLORset));
   memcpy(c2,c1,sizeof(COLORset));
   memcpy(c1,t,sizeof(COLORset));
}

static int COLORset_less(COLORset *c1,COLORset *c2)
{
   int i;
   if (c1->count != c2->count) {
      return c1->count < c2->count;
   }
   for (i = 0; i < c1->count;++i) {
      if (c1->members[i] != c2->members[i]) {
         return c1->members[i] < c2->members[i];
      }
   }
   return 0;
}

void COLORset_quicksort (COLORset *cclasses, int ccount)
{
    int i, j;
    COLORset temp,t;
    if (ccount <= 1) return;

    COLORset_SWAP (&(cclasses[0]), &(cclasses[(ccount - 1)/2]), &temp);

    i = 0; j = ccount;
    memcpy(&t,&(cclasses[0]),sizeof(COLORset));

    while (1) {
       do i++; while (i < ccount && COLORset_less(&(cclasses[i]),&t));
       do j--; while ( COLORset_less(&t,&(cclasses[j]) ) );
        if (j < i) break;
        COLORset_SWAP (& (cclasses[i]), &(cclasses[j]), &temp);
    }
    COLORset_SWAP (&(cclasses[0]), &(cclasses[j]), &temp);

    COLORset_quicksort (cclasses, j);
    COLORset_quicksort (cclasses + i, ccount - i);
}

static void COLORset_unify (COLORset *cclasses, int* new_ccount, int ccount)
{
   int i;
   COLORset temp;
   COLORset_quicksort(cclasses,ccount);

   *new_ccount = 0;
   i = 0;
   if (! ccount) return;

   /* Find first non-empty set */
   while(!cclasses[i].count) {
      COLORfree_set(&(cclasses[i++]));
   }

   for (; i < ccount; ++i) {
      if (*new_ccount == 0 || COLORset_less(&(cclasses[*new_ccount -1]),&(cclasses[i]))) {
         (*new_ccount)++;
         if(*new_ccount  < i + 1) {
            COLORset_SWAP(&(cclasses[*new_ccount -1]),& (cclasses[i]),&temp);
         }
      } else {
         COLORfree_set(&(cclasses[i]));
      }
   }
}


static int new_eindex(int v, int v1, int v2)
{
   if (v == v2) {return v1;}
   if (v > v2)  {return v - 1;}

   return v;
}

static int are_in_same(const COLORset* debugcolors,int ndebugcolors,int v1,int v2)
{
   int v1_color = -1;
   int v2_color = -2;
   int i;
   for (i = 0; i < ndebugcolors; ++i) {
      int j;
      for (j = 0; j < debugcolors[i].count; ++j) {
         if (debugcolors[i].members[j] == v1) {
            v1_color = i;
         }
         if (debugcolors[i].members[j] == v2) {
            v2_color = i;
         }
      }
   }
   return (v1_color == v2_color);
}

static int create_same (colordata* parent_cd, int v1, int v2)
{
   int rval = 0;
   int i;
   int v1_covered = 0;
   COLORadjgraph G;
   colordata*    cd = (colordata*) NULL;
   int* Nv1Nv2Intersect = (int*) COLOR_SAFE_MALLOC(parent_cd->ncount,int);
   COLORcheck_NULL(Nv1Nv2Intersect,"Failed to allocate Nv1Nv2Intersect");

   cd = (colordata*) COLOR_SAFE_MALLOC(1,colordata);
   COLORcheck_NULL(cd,"Failed to allocate cd");

   init_colordata(cd);
   cd->depth = parent_cd->depth + 1;
   parent_cd->same_child = cd;

   cd->v1 = v1;
   cd->v2 = v2;

   cd->upper_bound = parent_cd->upper_bound;
   cd->lower_bound = parent_cd->lower_bound;
   cd->dbl_lower_bound = parent_cd->dbl_lower_bound;


   /* Create Nv1Nv2Intersect. The intersection is labeled with >= 2
      (exactly two in a simple graph) */
   rval = COLORadjgraph_build(&G,parent_cd->ncount,parent_cd->ecount, parent_cd->elist);
   COLORcheck_rval(rval,"COLORadjgraph_build");

   for(i = 0; i < parent_cd->ncount; ++i) {
      Nv1Nv2Intersect[i] = 0;
   }
   for(i = 0; i < G.nodelist[v2].degree; ++i) {
      int v_i = G.nodelist[v2].adj[i];
      Nv1Nv2Intersect[v_i] = 1;
   }
   COLORadjgraph_free(&G);


   /* Create contracted graph */
   cd->ncount = parent_cd->ncount - 1;
   cd->ecount = parent_cd->ecount;
   cd->elist  = (int*) COLOR_SAFE_MALLOC(2 * cd->ecount,int);
   COLORcheck_NULL(cd->elist,"Failed to allocate cd->elist");
   memcpy(cd->elist,parent_cd->elist, 2 * parent_cd->ecount * sizeof(int));
   for(i = 0; i < cd->ecount; ++i) {
      cd->elist[2*i] = new_eindex(cd->elist[2*i],v1,v2);
      cd->elist[2*i+1] = new_eindex(cd->elist[2*i+1],v1,v2);
   }
   rval = COLORadjgraph_build(&G,cd->ncount,cd->ecount, cd->elist);
   COLORcheck_rval(rval,"COLORadjgraph_build");

   rval =  COLORadjgraph_simplify(&G);
   COLORcheck_rval(rval,"COLORadjgraph_simplify");

   rval = COLORadjgraph_extract_edgelist(&(cd->ecount), &(cd->elist),&G);

   cd->orig_node_ids = (int*) COLOR_SAFE_MALLOC(cd->ncount,int);
   COLORcheck_NULL(cd->orig_node_ids,"Failed to allocate cd->orig_node_ids");
   for (i = 0; i < cd->ncount; ++i) {
      int j = (i < v2 ) ?  i : i + 1;
      cd->orig_node_ids[i] = parent_cd->orig_node_ids[j];
   }
   cd->parent = parent_cd;
   cd->debugcolors = parent_cd->debugcolors;
   cd->ndebugcolors = parent_cd->ndebugcolors;

   /* END create contracted graph */

   if (COLORdbg_lvl() > 1) {
      printf("create_same created following graph:\n");
      COLORgraph_print(cd->ecount,cd->elist);
   }

   /* Transfer independent sets: */
   cd->gallocated = cd->ccount   =  parent_cd->ccount + 1;
   cd->cclasses = (COLORset*) COLOR_SAFE_MALLOC(cd->gallocated,COLORset);
   for (i = 0; i < parent_cd->ccount; ++i) {
      int j;
      int add_v1 = 1;
      cd->cclasses[i].members = (int*) COLOR_SAFE_MALLOC(parent_cd->cclasses[i].count,int);
      cd->cclasses[i].count = 0;
      for (j = 0; j < parent_cd->cclasses[i].count; ++j) {
         if (Nv1Nv2Intersect[parent_cd->cclasses[i].members[j]] == 1) {
            add_v1 = 0;
            j = parent_cd->cclasses[i].count;/*break*/
         }
      }
      for (j = 0; j < parent_cd->cclasses[i].count; ++j) {
         if (parent_cd->cclasses[i].members[j] == v1) {
            if (add_v1) {v1_covered = 1;}
            else {continue;}
         }
         if (parent_cd->cclasses[i].members[j] < v2) {
            cd->cclasses[i].members[(cd->cclasses[i].count)++] =
               parent_cd->cclasses[i].members[j];
         }
         else if (parent_cd->cclasses[i].members[j] >  v2) {
            cd->cclasses[i].members[(cd->cclasses[i].count)++] =
               parent_cd->cclasses[i].members[j] - 1;
         }
         /* else 'parent_cd->cclasses[i].members[j] == v2' and we skip it*/
      }
      if(COLORdbg_lvl() > 1) {
         printf("PARENT SET ");
         for (j = 0; j < parent_cd->cclasses[i].count; ++j) {
            printf(" %d",parent_cd->cclasses[i].members[j]);
         }
         printf("\n");
         printf("TRANS SET ");
         for (j = 0; j < cd->cclasses[i].count; ++j) {
            printf(" %d",cd->cclasses[i].members[j]);
         }
         printf("\n");
      }
   }
   if (!v1_covered) {
      /* Create the singular set v1 as last set, because we might not add
         v1 to any set: */
      printf("Adding extra set %d\n", v1);
      cd->cclasses[parent_cd->ccount].count  = 1;
      cd->cclasses[parent_cd->ccount].members = (int*) COLOR_SAFE_MALLOC(1,int);
      cd->cclasses[parent_cd->ccount].members[0] = v1;
   } else {
      cd->cclasses[parent_cd->ccount].count  = 0;
      cd->cclasses[parent_cd->ccount].members = (int*) 0;
      cd->ccount--;
   }
   /** prune duplicate  sets */
   {
      int j;
      COLORset_unify (cd->cclasses, &(cd->ccount), cd->ccount);
      for (i = 0 ; i < cd->ccount; ++i ) {
         if (COLORdbg_lvl() > 1) {
            printf("TRANSSORT SET ");
            for (j = 0; j < cd->cclasses[i].count; ++j) {
               printf(" %d",cd->cclasses[i].members[j]);
            }
            printf("\n");
         }
         rval = COLORcheck_set(&(cd->cclasses[i]),cd->ncount,cd->ecount,cd->elist);
         COLORcheck_rval(rval, "Illegal colorset created in create_same\n!");
      }
   }
   /* END Transfer independent sets: */

 CLEANUP:
   if (rval) {
      if (cd) {
         free_colordata(cd);
         free(cd);
      }
      parent_cd->same_child = (colordata*) NULL;
   }
   COLORadjgraph_free(&G);
   if (Nv1Nv2Intersect) free(Nv1Nv2Intersect);
   return rval;
}

static int create_differ(colordata* parent_cd, int v1, int v2)
{
   int rval = 0;
   int i;
   COLORadjgraph G;
   int           v2_covered;
   colordata*    cd = (colordata*) COLOR_SAFE_MALLOC(1,colordata);
   COLORcheck_NULL(cd,"Failed to allocate cd");


   init_colordata(cd);

   cd->depth = parent_cd->depth + 1;
   parent_cd->diff_child = cd;

   cd->v1 = v1;
   cd->v2 = v2;

   cd->upper_bound = parent_cd->upper_bound;
   cd->lower_bound = parent_cd->lower_bound;
   cd->dbl_lower_bound = parent_cd->dbl_lower_bound;


   /* Create  graph with extra edge (v1,v2) */
   cd->ncount = parent_cd->ncount;
   cd->ecount = parent_cd->ecount + 1;
   cd->elist  = (int*) COLOR_SAFE_MALLOC(2 * cd->ecount,int);
   COLORcheck_NULL(cd->elist,"Failed to allocate cd->elist");

   memcpy(cd->elist,parent_cd->elist, 2 * parent_cd->ecount * sizeof(int));
   cd->elist[ 2 * (cd->ecount - 1)] = v1;
   cd->elist[ 2 * (cd->ecount - 1) + 1] = v2;

   rval = COLORadjgraph_build(&G, cd->ncount,cd->ecount,cd->elist);
   COLORcheck_rval(rval,"COLORadjgraph_build failed");
   rval = COLORadjgraph_simplify(&G);
   COLORcheck_rval(rval,"COLORadjgraph_simplify failed");
   COLORadjgraph_extract_edgelist(&cd->ecount, &cd->elist,&G);
   COLORcheck_rval(rval,"COLORadjgraph_extract_edgelist");

   cd->orig_node_ids = (int*) COLOR_SAFE_MALLOC(cd->ncount,int);
   COLORcheck_NULL(cd->orig_node_ids,"Failed to allocate cd->orig_node_ids");
   for (i = 0; i < cd->ncount; ++i) {
      cd->orig_node_ids[i] = parent_cd->orig_node_ids[i];
   }
   cd->parent = parent_cd;
   cd->debugcolors = parent_cd->debugcolors;
   cd->ndebugcolors = parent_cd->ndebugcolors;

   /* END: Create  graph with extra edge (v1,v2) */

   if (COLORdbg_lvl() > 1) {
      printf("create_differ created following graph:\n");
      COLORgraph_print(cd->ecount,cd->elist);
   }

   /* Transfer independent sets by removing v2 if both v1 and v2 are currently contained: */
   cd->gallocated = cd->ccount   =  parent_cd->ccount + 1;
   cd->cclasses = (COLORset*) COLOR_SAFE_MALLOC(cd->gallocated,COLORset);

   for (i = 0; i < parent_cd->ccount; ++i) {
      int j;
      int v1_found = 0;
      cd->cclasses[i].members = (int*) COLOR_SAFE_MALLOC(parent_cd->cclasses[i].count,int);
      COLORcheck_NULL(cd->cclasses[i].members,"Failed to allocate cd->cclasses[i].members");
      cd->cclasses[i].count = 0;

      for (j = 0; j < parent_cd->cclasses[i].count; ++j) {
         int current_elm = parent_cd->cclasses[i].members[j];
         if (current_elm ==  v1) {
            v1_found = 1;
         }
         if (current_elm ==  v2) {
            if (v1_found) {
               continue;
            } else {
               v2_covered = 1;
            }
         }
         cd->cclasses[i].members[cd->cclasses[i].count] = current_elm;
         (cd->cclasses[i].count)++;
      }
      if(COLORdbg_lvl() > 1 ) {
         printf("PARENT SET ");
         for (j = 0; j < parent_cd->cclasses[i].count; ++j) {
            printf(" %d",parent_cd->cclasses[i].members[j]);
         }
         printf("\n");

         printf("TRANS SET ");
         for (j = 0; j < cd->cclasses[i].count; ++j) {
            printf(" %d",cd->cclasses[i].members[j]);
         }
         printf("\n");
      }
      rval = COLORcheck_set(&(cd->cclasses[i]),cd->ncount,cd->ecount,cd->elist);
      COLORcheck_rval(rval, "Illegal colorset created in create_same\n!");

   }
   if (!v2_covered) {
      /* Create the singular set v2 as last set, because we might not add
         v2 to any set: */
      cd->cclasses[parent_cd->ccount].count  = 1;
      cd->cclasses[parent_cd->ccount].members = (int*) COLOR_SAFE_MALLOC(1,int);
      cd->cclasses[parent_cd->ccount].members[0] = v2;
   } else {
      cd->cclasses[parent_cd->ccount].count   = 0;
      cd->cclasses[parent_cd->ccount].members = (int*) 0;
      cd->ccount--;
   }
   /* END Transfer independent sets: */


   /** prune duplicate  sets */
   {
      int j;
      COLORset_unify (cd->cclasses, &(cd->ccount), cd->ccount);
      for (i = 0 ; i < cd->ccount; ++i ) {
         if (COLORdbg_lvl() > 1) {
            printf("TRANSSORT SET ");
            for (j = 0; j < cd->cclasses[i].count; ++j) {
               printf(" %d",cd->cclasses[i].members[j]);
            }
            printf("\n");
         }
         rval = COLORcheck_set(&(cd->cclasses[i]),cd->ncount,cd->ecount,cd->elist);
         COLORcheck_rval(rval, "Illegal colorset created in create_same\n!");
      }
   }

 CLEANUP:
   if (rval) {
      if (cd) {
         free_colordata(cd);
         free(cd);
      }
      parent_cd->diff_child = (colordata*) NULL;
   }
   COLORadjgraph_free(&G);

   return rval;
}

static int grab_integral_solution(colordata* cd,
                                  double* x,
                                  double tolerance)
{
   int rval = 0;
   double incumbent;
   int* colored = (int*) NULL;
   int i;

   rval = COLORlp_objval (cd->lp, &incumbent);
   COLORcheck_rval (rval, "COLORlp_objval failed");

   colored = (int*) COLOR_SAFE_MALLOC(cd->ncount,int);
   COLORcheck_NULL(colored,"Failed to allocate colored");
   for (i = 0; i < cd->ncount; ++i) {colored[i] = 0;}

   COLORfree_sets(&(cd->bestcolors),&(cd->nbestcolors));
   cd->bestcolors = (COLORset*) realloc(cd->bestcolors,
                                        (int)incumbent * sizeof(COLORset));
   COLORcheck_NULL(cd->bestcolors,"Failed to realloc cd->bestcolors");

   cd->nbestcolors = 0;
   for (i = 0; i < cd->ccount; ++i) {
      if (x[i] >= 1.0  - tolerance) {
         int j = cd->nbestcolors;
         int k;
         cd->bestcolors[j].count = 0 ;
         cd->bestcolors[j].members = (int*) COLOR_SAFE_MALLOC(cd->cclasses[i].count,int);
         COLORcheck_NULL(cd->bestcolors[j].members,
                         "Failed to realloc cd->bestcolors[j].members");

         for(k = 0; k < cd->cclasses[i].count; ++k) {
            if (!colored[cd->cclasses[i].members[k]]) {
               colored[cd->cclasses[i].members[k]] = 1;
               cd->bestcolors[j].members[cd->bestcolors[j].count++] =
                  cd->cclasses[i].members[k];
            }
         }
         cd->nbestcolors++;
      }
   }
   rval = COLORcheck_coloring(cd->bestcolors, cd->nbestcolors,
                              cd->ncount, cd->ecount, cd->elist);
   COLORcheck_rval(rval,"ERROR: An incorrect coloring was created.");

   printf("Intermediate coloring:\n");
   print_colors(cd->bestcolors,cd->nbestcolors);
   assert(cd->nbestcolors == (int) incumbent);

   if (cd->nbestcolors < cd->upper_bound) {
      cd->upper_bound = cd->nbestcolors;
   }
   cd->status = finished;
 CLEANUP:
   if (colored) free(colored);
   return rval;
}

static int x_frac(const double x[], int i)
{
   double mean = 0.55;
   double frac = mean - fabs(x[i] - mean);
   assert(frac <= mean);
   return (int) ( frac * (double) INT_MAX);
}


/** compute array-index from row-index v1 and column-index v2.*/
static int nodepair_ref_key(int v1, int v2)
{
   assert(v1 < v2);/* We store only the elements of the lower left triangle within the ncount x ncount matrix */
   return v2 * (v2 - 1) / 2 + v1;
}

/** compute row-index v1 and column-index v2 from array-index.*/
static void inodepair_ref_key(int* v1, int* v2,int index)
{
   *v2 = (int) floor (sqrt(2*((double)index) + 0.25) + 0.5);
   *v1 = index - (*v2 * (*v2 -1) / 2);
}

static int insert_fractional_pairs_to_heap(colordata* cd, const double x[],
                                           int*          nodepair_refs,
                                           double*       nodepair_weights,
                                           int           npairs,
                                           COLORNWTHeap* heap)
{
   int rval = 0;
   int i;
   int ref_key;

   for (i = 0; i < cd->ccount; ++i) {
      int j;
      if (x[i] <= 0.0 || x[i] >= 1.0 ) {continue;}
      for (j = 0; j < cd->cclasses[i].count; ++j) {
         int v1 = cd->cclasses[i].members[j];
         int k;
         for (k = j + 1 ; k < cd->cclasses[i].count; ++k) {
            assert (k != j);
            int v2 = cd->cclasses[i].members[k];
            assert(v1 < v2);
            ref_key  = nodepair_ref_key(v1, v2);

            nodepair_weights[ref_key] += x[i];
            fflush(stdout);
         }
      }
   }

   for (ref_key = 0; ref_key < npairs; ++ref_key) {
      if (nodepair_weights[ref_key] > 0.0) {
         int         heap_key  =  - x_frac(nodepair_weights,ref_key);

         rval = COLORNWTheap_insert(heap,
                                    &(nodepair_refs[ref_key]),
                                    heap_key,
                                    & (nodepair_refs[ref_key]));
         COLORcheck_rval(rval, "Failed in COLORNWTheap_insert");
      }
   }

 CLEANUP:
   return rval;
}


static int find_strongest_children(int           *strongest_v1,
                                   int           *strongest_v2,
                                   colordata*    cd,
                                   COLORNWTHeap* cand_heap,
                                   int*          nodepair_refs,
                                   double*       nodepair_weights)
{
   int    rval = 0;

   int    max_non_improving_branches  = cd->ncount / 100 + 1;
   int    remaining_branches          = max_non_improving_branches;
   double strongest_dbl_lb = 0.0;
   int*   min_nodepair;
   *strongest_v1 = -1;
   *strongest_v2 = -1;

   while ( (min_nodepair = (int*) COLORNWTheap_min(cand_heap)) && (remaining_branches--) ) {
      int v1 = -1,v2 = -1;
      double dbl_child_lb = (double) cd->ncount;
      inodepair_ref_key(&v1,&v2, (int)(min_nodepair - nodepair_refs));

      assert(v1 < v2);
      printf("Creating branches for v1 = %d, v2 = %d (node-pair weight %f)\n",v1,v2,
             nodepair_weights[(int)(min_nodepair-nodepair_refs)]);
      /* Create DIFFER and SAME */
      rval = create_same(cd,v1,v2);
      COLORcheck_rval(rval, "Failed in create_same");

      rval = create_differ(cd,v1,v2);
      COLORcheck_rval(rval, "Failed in create_differ");

      compute_lower_bound(cd->same_child);
      compute_lower_bound(cd->diff_child);

      dbl_child_lb = (cd->same_child->dbl_lower_bound < cd->diff_child->dbl_lower_bound) ?
         cd->same_child->dbl_lower_bound : cd->diff_child->dbl_lower_bound;

      free_colordata(cd->same_child);
      free(cd->same_child); cd->same_child = (colordata*) NULL;
      free_colordata(cd->diff_child);
      free(cd->diff_child); cd->diff_child = (colordata*) NULL;

      if (dbl_child_lb > strongest_dbl_lb) {
         strongest_dbl_lb = dbl_child_lb;
         *strongest_v1     = v1;
         *strongest_v2     = v2;
         remaining_branches = max_non_improving_branches;
      }
      printf("Found child bound of %f for v1 = %d, v2 = %d, nodepair_weight = %f .\n",dbl_child_lb,
             v1, v2,
             nodepair_weights[(int)(min_nodepair-nodepair_refs)]);
   }
   {
      int nodepair_ref = nodepair_ref_key(*strongest_v1,*strongest_v2);
      printf("Found strongest child bound of %f for v1 = %d, "
             "v2 = %d, nodepair_weight = %f .\n",
             strongest_dbl_lb, *strongest_v1,
             *strongest_v2, nodepair_weights[nodepair_ref]);
   }
 CLEANUP:
   return rval;
}

static int create_branches(colordata* cd)
{
   int rval = 0;
   int i;
   double* x = (double*) NULL;
   int  strongest_v1 = -1, strongest_v2 = -1;
   /* For collecting node pairs we simply mark collected pairs in a
      upper triangle matrix.  The smaller id determines the row and
      the larger the column.
   */
   int*    nodepair_refs    = (int*) NULL;
   double* nodepair_weights = (double*) NULL;
   int     npairs = cd->ncount * (cd->ncount - 1) / 2;
   /* For each vertex marked in the nodepair bucket
      we mark its most fractional column in s1_value.*/
   int*    mf_col = (int*) NULL;
   COLORNWTHeap*  cand_heap = (COLORNWTHeap*) NULL;

   printf("Entered create_branches\n");

   rval = COLORNWTheap_init(&cand_heap,npairs);
   COLORcheck_rval(rval,"Failed in COLORNWTheap_init");

   nodepair_refs = (int*) COLOR_SAFE_MALLOC(npairs, int);
   COLORcheck_NULL(nodepair_refs,"Failed ot allocate nodepair_refs");
   for (i = 0; i < npairs; ++i) {nodepair_refs[i] = -1;}


   nodepair_weights = (double*) COLOR_SAFE_MALLOC(npairs, double);
   COLORcheck_NULL(nodepair_weights,"Failed ot allocate nodepair_weights");
   for (i = 0; i < npairs; ++i) {nodepair_weights[i] = 0;}


   mf_col = (int*) COLOR_SAFE_MALLOC(cd->ncount, int);
   COLORcheck_NULL(mf_col,"Failed ot allocate nodepair_bucket");
   for (i = 0; i < cd->ncount; ++i) {mf_col[i] = -1;}

   x = (double*) COLOR_SAFE_MALLOC(cd->ccount,double);
   COLORcheck_NULL(x,"Failed ot allocate x");

   rval = COLORlp_optimize(cd->lp);
   COLORcheck_rval (rval, "COLORlp_optimize failed");

   rval = COLORlp_x(cd->lp,x);
   COLORcheck_rval(rval,"Failed in COLORlp_x");

   rval = insert_fractional_pairs_to_heap(cd, x,nodepair_refs,
                                          nodepair_weights,npairs,
                                          cand_heap);
   COLORcheck_rval(rval, "Failed in insert_fractional_paris_to_heap");

   if (COLORNWTheap_size(cand_heap) == 0) {
      printf("LP returned integral solution.\n");
      grab_integral_solution(cd,x,0.0);
      goto CLEANUP;
   }

   printf("Collected %d branching candidates.\n",COLORNWTheap_size(cand_heap));

   rval = find_strongest_children(&strongest_v1,&strongest_v2,
                                  cd,cand_heap,nodepair_refs,
                                  nodepair_weights);
   COLORcheck_rval(rval, "Failed in find_strongest_children");

   /* Create DIFFER and SAME for strongest children */
   rval = create_same(cd,strongest_v1,strongest_v2);
   COLORcheck_rval(rval, "Failed in create_same");

   rval = create_differ(cd,strongest_v1,strongest_v2);
   COLORcheck_rval(rval, "Failed in create_differ");





   if (cd->same_child && cd->ndebugcolors) {
      int same_opt_track = 0;

      if (cd->opt_track) {
         same_opt_track = are_in_same(cd->debugcolors,cd->ndebugcolors,
                                      cd->orig_node_ids[cd->same_child->v1],
                                      cd->orig_node_ids[cd->same_child->v2]);
         cd->same_child->opt_track = same_opt_track;
         cd->diff_child->opt_track = !same_opt_track;
      } else  {
         cd->same_child->opt_track = 0;
         cd->diff_child->opt_track = 0;
      }
   }

 CLEANUP:
   free_lbcolordata(cd);

   if (cand_heap) {
      COLORNWTheap_free(cand_heap);
      cand_heap = (COLORNWTHeap*) NULL;
   }
   if(x)                free (x);
   if(mf_col)           free (mf_col);
   if(nodepair_refs)    free (nodepair_refs);
   if(nodepair_weights) free (nodepair_weights);

   return rval;
}

/** Transform the coloring of cd->same_child into a coloring of cd.
 */
static int collect_same_child(colordata* cd)
{
   int rval = 0;

   if (cd->same_child && cd->same_child->nbestcolors) {
      if (!cd->nbestcolors || cd->same_child->nbestcolors < cd->upper_bound) {
         int i;
         if (cd->nbestcolors) {
            COLORfree_sets(&(cd->bestcolors),&(cd->nbestcolors));
         }
         cd->upper_bound = cd->nbestcolors = cd->same_child->nbestcolors;
         cd->same_child->nbestcolors = 0;
         cd->bestcolors = cd->same_child->bestcolors;
         cd->same_child->bestcolors = (COLORset*) NULL;
         for (i = 0; i < cd->nbestcolors; ++i) {
            int j;
            int add_v2 = 0;
            for (j = 0; j < cd->bestcolors[i].count ; ++j) {
               if (cd->bestcolors[i].members[j] == cd->same_child->v1) {
                  assert(add_v2 == 0); add_v2 = 1;
               }
               if (cd->bestcolors[i].members[j] >= cd->same_child->v2) {
                  (cd->bestcolors[i].members[j])++;
               }
            }
            if (add_v2) {
               int k ;
               (cd->bestcolors[i].count)++;
               cd->bestcolors[i].members =
                  (int*) realloc(cd->bestcolors[i].members,
                                 cd->bestcolors[i].count * sizeof(int));
               COLORcheck_NULL(cd->bestcolors[i].members,
                               "Failed to realloc cd->bestcolors[i].members");
               k = cd->bestcolors[i].count - 1;
               while ( (k > 0) &&  cd->bestcolors[i].count &&
                       (cd->same_child->v2 < cd->bestcolors[i].members[k - 1]) ) {
                  k--;
               }
               if (k < cd->bestcolors[i].count - 1) {
                  memmove(cd->bestcolors[i].members + k + 1,
                          cd->bestcolors[i].members + k,
                          (cd->bestcolors[i].count - 1 - k) * sizeof(int));
               }
               cd->bestcolors[i].members[k] = cd->same_child->v2;

            }
         }
         print_colors(cd->bestcolors,cd->nbestcolors);
         fflush(stdout);
         rval = COLORcheck_coloring(cd->bestcolors,cd->nbestcolors,
                                    cd->ncount, cd->ecount, cd->elist);
         COLORcheck_rval(rval,"ERROR: An incorrect coloring was created.");
      }
   }
 CLEANUP:
   return rval;
}

/** Transform the coloring of cd->same_child into a coloring of cd.
 */
static int collect_diff_child(colordata* cd)
{
   int rval = 0;

   if (cd->diff_child && cd->diff_child->nbestcolors) {
      if (!cd->nbestcolors || cd->diff_child->nbestcolors < cd->upper_bound) {
         if (cd->nbestcolors) {
            COLORfree_sets(&(cd->bestcolors),&(cd->nbestcolors));
         }
         cd->upper_bound = cd->nbestcolors = cd->diff_child->nbestcolors;
         cd->diff_child->nbestcolors = 0;
         cd->bestcolors = cd->diff_child->bestcolors;
         cd->diff_child->bestcolors = (COLORset*) NULL;

         rval = COLORcheck_coloring(cd->bestcolors,cd->nbestcolors,
                                    cd->ncount, cd->ecount, cd->elist);
         COLORcheck_rval(rval,"ERROR: An incorrect coloring was created.");
      }
   }

 CLEANUP:
   return rval;
}

static int print_graph_operations(const colordata* cd)
{
   if (cd->parent) {
      print_graph_operations(cd->parent);
      if (cd->ncount < cd->parent->ncount)
         printf("SAME ");
      else
         printf("DIFF ");
      printf("%d %d\n",
             cd->parent->orig_node_ids[cd->v1],
             cd->parent->orig_node_ids[cd->v2]);
   }
   return 0;
}

static int compute_lower_bound(colordata* cd)
{
   int rval = 0;
   int iterations       = 0;
   int break_while_loop = 1;

   double last_snapshot_time;

   double cur_time;
   double lb_rtime;
/*    int    parent_lb = cd->lower_bound; */

   printf("Starting compute_lower_bound with lb %d (%f) and ub %d at depth %d (id = %d, opt_track = %d).\n",
          cd->lower_bound,cd->dbl_lower_bound,cd->upper_bound,cd->depth,cd->id, cd->opt_track);


   lb_rtime = COLORcpu_time();
   last_snapshot_time = lb_rtime;

   reset_ages(cd->cclasses,cd->ccount) ;
   cd->gallocated = cd->ccount;

   rval = build_lp(cd);
   COLORcheck_rval (rval, "build_lp failed");

   rval = COLORstable_initenv (&(cd->mwis_env),cd->pname,write_mwis);
   COLORcheck_rval (rval, "COLORgreedy failed");


   cd->mwis_pi = (COLORNWT *) COLOR_SAFE_MALLOC (cd->ncount,COLORNWT);
   COLORcheck_NULL (cd->mwis_pi, "out of memory for mwis_pi");

   cd->retirementage = cd->ncount / 4 + 1;
   do {
      ++iterations;

      if (iterations > cd->retirementage && cd->ccount > 3 * cd->ncount) {
         delete_old_colorclasses(cd);
      }

      if (cclasses_outfile) {
         int add_timestamp = 1;
         cur_time = COLORcpu_time();
         if (cur_time - last_snapshot_time > 7200) {
            write_snapshot(cd,add_timestamp);
            last_snapshot_time = cur_time;
         }
      }
      rval = COLORlp_optimize(cd->lp);
      COLORcheck_rval (rval, "COLORlp_optimize failed");

      rval = grow_ages(cd);
      COLORcheck_rval (rval, "grow_ages failed");

      if (COLORdbg_lvl() > 1) {print_ages(cd);}

      rval = COLORlp_pi (cd->lp, cd->pi);
      COLORcheck_rval (rval, "COLORlp_pi failed");

      make_pi_feasible(cd);

      COLOR_double2COLORNWT(cd->mwis_pi,&(cd->mwis_pi_scalef),
                            cd->pi,cd->ncount);

      compute_objective(cd);

/*       if ((double) parent_lb < cd->dbl_lower_bound) { */
      {
         int set_i;

         rval = COLORstable_wrapper(&(cd->mwis_env),&(cd->newsets), &(cd->nnewsets),
                                    cd->ncount, cd->ecount,
                                    cd->elist, cd->mwis_pi,cd->mwis_pi_scalef);
         COLORcheck_rval (rval, "COLORstable_gurobi failed");

         for (set_i = 0; set_i < cd->nnewsets; ++set_i) {
            rval = COLORlp_addcol (cd->lp, cd->newsets[set_i].count,
                                   cd->newsets[set_i].members,
                                   cd->coef, 1.0, 0.0, 1.0,
                                   COLORlp_CONTINUOUS, NULL);
         }
         break_while_loop = (cd->nnewsets == 0);

         concat_newsets(cd);
      }

   } while ( (iterations < cd->maxiterations) &&
/*              ((double) parent_lb < cd->dbl_lower_bound) && */
             !break_while_loop);
   if (iterations < cd->maxiterations) {

      printf ("Found bound of %lld (%20.16g), upper_bound %d (id = %d, iterations = %d, opt_track = %d).\n",
              (long long) cd->lower_bound,cd->dbl_lower_bound,
              cd->upper_bound,cd->id,iterations, cd->opt_track);

      if (COLORdbg_lvl()) {
         print_graph_operations(cd);
      }

      delete_old_colorclasses(cd);

      lb_rtime = COLORcpu_time() - lb_rtime;
      cd->status = LP_bound_computed;
      printf("Computing initial lower bound took %f seconds.\n",lb_rtime);
   } else {
      fprintf (stderr, "ERROR: Lower bound could not be found in %d iterations!\n", iterations);
      rval = 1; goto CLEANUP;
   }
   fflush(stdout);

 CLEANUP:
   return rval;
}





static int insert_into_branching_heap(COLORNWTHeap* heap, colordata* cd, double key_mult)
{
   int rval = compute_lower_bound (cd);
   int dummy_href;
   /** We use the lower bound adjusted by the depth of a node as a
       heap key.  The idea behind this is: in case of equal bounds
       deeper nodes should be preferred to potentially improve the
       upper bound faster.
    */
   COLORNWT heap_key = (COLORNWT) (cd->dbl_lower_bound * key_mult) - cd->depth;
   COLORcheck_rval(rval,"Failed in compute_lower_bound (cd);");
   rval = COLORNWTheap_insert(heap,&dummy_href,
                              heap_key,
                              (void*) cd);
   COLORcheck_rval(rval, "Failed to COLORNWTheap_insert");
 CLEANUP:
   return rval;
}

static int remove_finished_subtree(colordata* child)
{
   int rval = 0;
   colordata* cd = (colordata*) child;

   while (cd) {
      if (cd->same_child && cd->same_child->status == finished) {
         rval = collect_same_child(cd);
         COLORcheck_rval(rval,"Failed in collect_same_child");
         free_colordata(cd->same_child);
         free(cd->same_child);
         cd->same_child = (colordata*) NULL;
      }

      if (cd->diff_child && cd->diff_child->status == finished) {
         rval = collect_diff_child(cd);
         COLORcheck_rval(rval,"Failed in collect_diff_child");
         free_colordata(cd->diff_child);
         free(cd->diff_child);
         cd->diff_child = (colordata*) NULL;
      }
      if (!cd->same_child && !cd->diff_child) {
         cd->status = finished;
         cd = cd->parent;
      } else {
         cd = (colordata*) NULL;
      }
   }
 CLEANUP:
   return rval;
}

static int compute_coloring(colordata* root_cd)
{
   int           rval = 0;
   COLORNWTHeap* br_heap            = (COLORNWTHeap*) NULL;
   colordata*    cd;
   double        key_mult           = (double) (COLORNWT_MAX - 1) / root_cd->ncount;
   int           global_upper_bound = INT_MAX;

   COLORNWTheap_init(&br_heap,root_cd->ncount * root_cd->ncount);

   rval = insert_into_branching_heap(br_heap,root_cd,key_mult);
   COLORcheck_rval(rval,"Failed in insert_into_branching_heap");

   rval = write_mwis_instances(root_cd);
   COLORcheck_rval(rval,"Failed ing write_mwis_instances");

   if (cclasses_outfile) {
      int add_timestamp = 0;
      write_snapshot(root_cd,add_timestamp);
      COLORcheck_rval(rval,"Failed ing write_snapshot");
   }

   double colheur_rtime = COLORcpu_time();
   heur_colors_with_stable_sets(root_cd);
   colheur_rtime = COLORcpu_time() - colheur_rtime;

   global_upper_bound = (global_upper_bound > root_cd->upper_bound) ?
      root_cd->upper_bound :global_upper_bound;

   printf("Upper bound heuristic on root node took %f seconds.\n",colheur_rtime);



   while ( (cd = (colordata*) COLORNWTheap_min(br_heap) ) ) {
      cd->upper_bound = global_upper_bound;
      if (cd->lower_bound < cd->upper_bound) {
         printf("Branching with lb %d (%f) and ub %d at depth %d (id = %d, "
                "opt_track = %d, unprocessed nodes = %d).\n",
                cd->lower_bound,cd->dbl_lower_bound,cd->upper_bound,
                cd->depth,
                cd->id, cd->opt_track,COLORNWTheap_size(br_heap) );

         /* Create children. If the current LP-solution turns out to be intergal,
            cd->upper_bound might decrease!
          */
         rval = create_branches(cd);
         COLORcheck_rval(rval,"Failed to create_branches");


         if (cd->same_child) {
            printf("Adding same child to heap\n");
            rval = insert_into_branching_heap(br_heap,cd->same_child,key_mult);
            COLORcheck_rval(rval,"Failed in insert_into_cb_heap");
         }

         if (cd->diff_child) {
            printf("Adding diff child to heap\n");
            rval = insert_into_branching_heap(br_heap,cd->diff_child,key_mult);
            COLORcheck_rval(rval,"Failed in insert_into_cb_heap");
         }
      }
      assert (cd->lower_bound <= cd->upper_bound);

      if (global_upper_bound > cd->upper_bound) {
         global_upper_bound = cd->upper_bound;
      }
      
      /** This is not simply an else-branch because the cd->upper_bound
          can be decreased in create_branches if the current LP-relaxation is
          intergal.
      */
      if (cd->lower_bound == cd->upper_bound) {
         remove_finished_subtree(cd);
      }
   }

 CLEANUP:
   COLORNWTheap_free(br_heap); br_heap = (COLORNWTHeap*) NULL;

   return rval;
}

int main (int ac, char **av)
{
    int rval = 0;
    int i;
    double start_time;
    double tot_rtime;

    colordata  _colordata;
    colordata* cd = &_colordata;


    COLORset*  debugcolors = (COLORset*) NULL;
    int       ndebugcolors = 0;


    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    init_colordata(cd);


    get_problem_name(cd->pname,edgefile);


    if (debug) printf ("Debugging turned on\n");
    if (outfile) printf ("Output File: %s\n", outfile);
    fflush (stdout);


    rval = COLORread_dimacs (edgefile, &(cd->ncount), &(cd->ecount),
                             &(cd->elist), (int **) NULL);
    COLORcheck_rval (rval, "COLORread_diamcs failed");
    cd->orig_node_ids = (int*) COLOR_SAFE_MALLOC(cd->ncount,int);
    COLORcheck_NULL(cd->orig_node_ids,"Failed to allocate cd->orig_node_ids");
    for (i = 0; i < cd->ncount; ++i) {cd->orig_node_ids[i] = i;}
    start_time = COLORcpu_time();

    if (cclasses_infile != (char*) NULL) {
       rval = COLORstable_read_stable_sets(&(cd->cclasses),&(cd->ccount),
                                           cd->ncount,cclasses_infile,cd->pname);
       COLORcheck_rval(rval,"Failed in COLORstable_read_stable_sets");
    } else {
       rval = COLORgreedy (cd->ncount, cd->ecount, cd->elist,
                           &(cd->ccount), &(cd->cclasses));
       COLORcheck_rval (rval, "COLORgreedycd failed");

       /*     rval = COLORplot_graphviz(ncount,ecount,elist,0); */
       printf ("Greedy Colors: %d\n", cd->ccount); fflush (stdout);
       print_colors(cd->cclasses,cd->ccount);
       COLORcopy_sets(&(cd->bestcolors),&(cd->nbestcolors),
                      cd->cclasses,cd->ccount);
       cd->upper_bound = cd->nbestcolors < cd->upper_bound ? cd->nbestcolors : cd->upper_bound;
    }

    if (color_infile != (char*) NULL) {
       rval = COLORstable_read_stable_sets(&debugcolors,&ndebugcolors,
                                           cd->ncount,color_infile,cd->pname);
       COLORcheck_rval(rval,"Failed in COLORstable_read_stable_sets");
       rval = COLORcheck_coloring(debugcolors,ndebugcolors,cd->ncount,
                           cd->ecount,cd->elist);
       COLORcheck_rval(rval,"Failed in COLORcheck_coloring");
       cd->debugcolors = debugcolors;
       cd->ndebugcolors = ndebugcolors;
       cd->opt_track = 1;
    }


    rval = compute_coloring(cd);
    COLORcheck_rval(rval, "Failed to compute_coloring");


    printf ("Opt Colors: %d\n", cd->nbestcolors); fflush (stdout);
    print_colors(cd->bestcolors,cd->nbestcolors);


    tot_rtime = COLORcpu_time() - start_time;
    printf("Computing coloring took %f seconds.\n",tot_rtime);



CLEANUP:
    if (debugcolors) free (debugcolors);
    free_colordata(cd);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int rval = 0;

    while ((c = getopt (ac, av, "dmo:r:w:c:u:")) != EOF) {
        switch (c) {
        case 'd':
           /* each -d increases the verbosity by one.*/
           ++debug;
            break;
        case 'o':
            outfile = optarg;
            break;
        case 'm':
           write_mwis = 1;
           break;
        case 'r':
           cclasses_infile = optarg;
           break;
        case 'w':
           cclasses_outfile = optarg;
           break;
        case 'c':
           color_infile = optarg;
           break;
        case 'u':
           initial_upper_bound = atoi(optarg);
           break;
        default:
            usage (av[0]);
            rval = 1;  goto CLEANUP;
        }
    }

    if (ac <= optind) {
        rval = 1; goto CLEANUP;
    } else {
        edgefile = av[optind++];
    }

CLEANUP:

    if (rval) usage (av[0]);
    return rval;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage %s: [-see below-] edge_file\n", f);
    fprintf (stderr, "   -d    turn on debugging\n");
    fprintf (stderr, "   -o f  write coloring to file f\n");
    fprintf (stderr, "   -m    write final stable set and clique instances\n");
    fprintf (stderr, "   -r f  read initial stable sets from file f\n");
    fprintf (stderr, "   -w f  write stable sets from file f\n");
    fprintf (stderr, "   -c f  read initial coloring from file f\n");
    fprintf (stderr, "   -u int  initial upper bound f\n");
}

 static int get_problem_name(char* pname,const char* efname)
{
   int    rval = 0;
   int    len = 0;
   const char * fname = strrchr(efname,'/');
   char * lastdot = strrchr(efname,'.');
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
      rval = 1;
   }
   printf("Extracted problem name %s\n",pname);

   return 0;
}

double COLORwall_time (void)
{
    return (double) time (0);
}

double COLORcpu_time (void)
{
    struct rusage ru;
    double t;

    getrusage (RUSAGE_SELF, &ru);

    t = ((double) ru.ru_utime.tv_sec) +
        ((double) ru.ru_utime.tv_usec) / 1000000.0;
    return t;
}

