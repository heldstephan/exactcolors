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

#include "lp.h"
#include "graph.h"
#include "color.h"

#include "mwis.h"
#include "plotting.h"

static char *edgefile = (char *) NULL;
static char *outfile = (char *) NULL;
static char *cclasses_outfile = (char *) NULL;
static char *cclasses_infile = (char *) NULL;


static int debug = 0;
static int write_mwis = 0;
typedef struct colordata colordata;
struct colordata {
   /* The instance graph */
   int ncount;
   int ecount;
   int *elist;


   /* The column generation LP. */
   COLORlp * lp;
   double *coef;
   double *pi;

   COLORNWT  mwis_pi_scalef;
   COLORNWT *mwis_pi;
   COLORNWT  lower_bound;
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

   int maxiterations;
   int retirementage;

   int v1,v2;
   colordata* same_child;
   colordata* diff_child;
   
   char  pname[256];
};



int main (int ac, char **av);
static int parseargs (int ac, char **av);
static void usage (char *f);
static int get_problem_name(char* pname,const char* efname);
static int build_lp(colordata* cd);
static int compute_coloring(colordata* cd, int lb, int ub, int depth);
static int grab_integral_solution(colordata* cd, double* x);

static void init_colordata(colordata* cd)
{
   sprintf(cd->pname,"colorprob");
   cd->ncount = 0;
   cd->ecount = 0;
   cd->elist = (int *) NULL;

   cd->coef = (double *) NULL;
   cd->pi = (double *) NULL;
   cd->mwis_pi_scalef = 1.0;
   cd->mwis_pi = (COLORNWT *) NULL;
   cd->lower_bound = 0;
   cd->lower_scaled_bound = 0;
    
   cd->lp       = (COLORlp *) NULL;
   cd->mwis_env  = (MWISenv*) NULL;

   cd->newsets = (COLORset*) NULL;   
   cd->nnewsets = 0;

   cd->ccount = 0;
   cd->cclasses = (COLORset *) NULL;
   cd->gallocated = 0;     

   cd->bestcolors = (COLORset *) NULL;
   cd->nbestcolors = 0;

   cd->maxiterations = 1000000;
   cd->retirementage = 1000000;

   cd->v1 = cd->v2 = -1;
   cd->same_child = (colordata*)NULL;
   cd->diff_child = (colordata*)NULL;
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

   if (cd->lp) COLORlp_free (&(cd->lp));
   if (cd->elist) free (cd->elist);
   if (cd->coef) free (cd->coef);
   if (cd->pi) free (cd->pi);
   if (cd->mwis_pi) free(cd->mwis_pi);


    COLORfree_sets(&(cd->newsets),&(cd->nnewsets));
    COLORfree_sets(&(cd->cclasses),&(cd->gallocated));
    COLORfree_sets(&(cd->bestcolors),&(cd->nbestcolors));
    COLORstable_freeenv(&(cd->mwis_env));

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

static int write_snapshot(colordata* cd, int add_timestamp) 
{
   int rval = 0;
   if (cclasses_outfile != (char*) NULL) {
      char   fname[256];

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
                                           fname,cd->pname); 
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

    rval = COLORlp_write (cd->lp, "look.lp");
    COLORcheck_rval (rval, "COLORlp_write failed");

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

MAYBE_UNUSED static int heur_colors_with_stable_sets(colordata* cd)
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

   rval = COLORlp_optimize(cd->lp);
   COLORcheck_rval (rval, "COLORlp_optimize failed");


   rval =  COLORlp_x (cd->lp, colsol);
   COLORcheck_rval (rval, "COLORlp_x failed");

   rval = COLORlp_objval (cd->lp, &incumbent);
   COLORcheck_rval (rval, "COLORlp_objval failed");

   grab_integral_solution(cd,colsol);

   COLORlp_set_all_coltypes(cd->lp,GRB_CONTINUOUS);
   COLORcheck_rval (rval, "COLORlp_set_all_coltypes");

   printf ("Found lower bound of %lld and upper bound of %g. Current best coloring:\n", 
           (long long) cd->lower_bound, incumbent);

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

static int create_same(colordata* parent_cd, int v1, int v2)
{
   int rval = 0;
   int i;
   int v1_covered = 0;
   graph G;
   colordata* cd = (colordata*) NULL;
   int* Nv1Nv2Interset = (int*) COLOR_SAFE_MALLOC(parent_cd->ncount,int);
   COLORcheck_NULL(Nv1Nv2Interset,"Failed to allocate Nv1Nv2Interset");

   cd = (colordata*) COLOR_SAFE_MALLOC(1,colordata);
   COLORcheck_NULL(cd,"Failed to allocate cd");

   init_colordata(cd);
   parent_cd->same_child = cd;

   cd->v1 = v1;
   cd->v2 = v2;

   /* Create Nv1Nv2Interset. The intersection is labeled with >= 2
      (exactly two in a simple graph) */
   rval = COLORadjgraph_build(&G,parent_cd->ncount,parent_cd->ecount, parent_cd->elist);
   COLORcheck_rval(rval,"COLORadjgraph_build");
   
   for(i = 0; i < parent_cd->ncount; ++i) {
      Nv1Nv2Interset[i] = 0;
   }
   for(i = 0; i < G.nodelist[v2].degree; ++i) {
      int v_i = G.nodelist[v2].adj[i];
      Nv1Nv2Interset[v_i] = 1;
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
   /* END create contracted graph */

   if (COLORdbg_lvl() > 2) {
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
         if (Nv1Nv2Interset[parent_cd->cclasses[i].members[j]] == 1) {
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
      if(COLORdbg_lvl() ) {
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
         printf("TRANSSORT SET ");
         for (j = 0; j < cd->cclasses[i].count; ++j) {
            printf(" %d",cd->cclasses[i].members[j]);
         }
         printf("\n");
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
   if (Nv1Nv2Interset) free(Nv1Nv2Interset);
   return rval;
}

static int create_differ(colordata* parent_cd, int v1, int v2)
{
   int rval = 0;
   int i;
   graph G;
   int v2_covered;
   colordata* cd = (colordata*) COLOR_SAFE_MALLOC(1,colordata);
   COLORcheck_NULL(cd,"Failed to allocate cd");

   init_colordata(cd);
   parent_cd->diff_child = cd;

   cd->v1 = v1;
   cd->v2 = v2;

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
    
   /* END: Create  graph with extra edge (v1,v2) */

   if (COLORdbg_lvl() > 2) {
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
      if(COLORdbg_lvl() ) {
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
         printf("TRANSSORT SET ");
         for (j = 0; j < cd->cclasses[i].count; ++j) {
            printf(" %d",cd->cclasses[i].members[j]);
         }
         printf("\n");
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

static int grab_integral_solution(colordata* cd, double* x)
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
      if (x[i] > 1.0  - DBL_EPSILON) {
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
   printf("Intermediate coloring:\n");
   print_colors(cd->bestcolors,cd->nbestcolors);
   assert(cd->nbestcolors == (int) incumbent);
 CLEANUP:
   if (colored) free(colored);
   return rval;
}

static int create_branches(colordata* cd)
{
   int rval = 0;
   int i;
   int v1,v2 = -1;  /* The branching nodes*/
   double* x;
   double most_frac_x = 0.5;
   int    s1 = -1, s2 = -1; /* The branching columns.*/

   printf("Entered create_branches\n");

   x = (double*) COLOR_SAFE_MALLOC(cd->ccount,double);
   COLORcheck_NULL(x,"Failed ot allocate x");

   COLORlp_write(cd->lp,"crbra_debug.lp");
   
   rval = COLORlp_optimize(cd->lp);
   COLORcheck_rval (rval, "COLORlp_optimize failed");

   rval = COLORlp_x(cd->lp,x);
   COLORcheck_rval(rval,"Failed in COLORlp_x");
   
   /** Find most fractional column s1: */
   for (i = 0; i < cd->ccount; ++i) {
      double frac = fabs(x[i] - 0.5);
      if ( frac < most_frac_x) {
         most_frac_x = frac;
         s1 = i;
      }
      if (most_frac_x == 0.0) {
         break;
      }
   }
   if (s1 == -1) {
      printf("LP returned integral solution.\n");
      grab_integral_solution(cd,x);
      goto CLEANUP;
   }
   v1 = cd->cclasses[s1].members[0];

   /** Find other column, containing v1. Note: such a column must
       exist as s1 was fractional (0 < x_{s1} < 1) and v1 must be
       covered by columns of value >= 1.
   */
   for (i = 0; i < cd->ccount; ++i) {
      int j;
      if (i == s1) {continue;}
      /* Todo replace by binary search:*/
      for (j = 0 ; j < cd->cclasses[i].count; ++j) {
         if (cd->cclasses[i].members[j] > v1) {
            break;
         }
         if (cd->cclasses[i].members[j] == v1) {
            s2 = i;
            i  = cd->ccount; 
            break;
         }
      }
   }   
   assert(s2 >=0);
   /** Find v2:*/
   for (i = 0; 
        i < cd->cclasses[s1].count && 
           i < cd->cclasses[s2].count; 
        ++i) {
      if (cd->cclasses[s1].members[i] != cd->cclasses[s2].members[i]) {         
         if (cd->cclasses[s1].members[i] != v1) {
            v2 = cd->cclasses[s1].members[i];
         }
         else if (v2 == -1) {
            v2 = cd->cclasses[s2].members[i];
         } else {
            assert(0);
         }
         break;
      }
   }
   if (v2 == -1) {/* if v2 wasn't found set to larger element of larger set:*/
      printf("Setting to larger elm\n");fflush(stdout);
      v2 = cd->cclasses[s1].members[cd->cclasses[s1].count - 1];
      if (cd->cclasses[s1].count < cd->cclasses[s2].count) {
         v2 = cd->cclasses[s2].members[cd->cclasses[s2].count - 1];
      }
   }
   if (v2 < v1) {
      int t = v1; v1 = v2; v2 = t;
   }
   assert(v1 != v2);
   printf("Creating branches for v1 = %d, v2 = %d\n",v1,v2);
   /* Create DIFFER and SAME*/
   create_same(cd,v1,v2);
   COLORcheck_rval(rval, "Failed in create_same");

   rval = create_differ(cd,v1,v2);
   COLORcheck_rval(rval, "Failed in create_differ");

 CLEANUP:
   if(x) free(x);

   return rval;
}


static int collect_same_child(colordata* cd)
{
   int rval = 0;

   if (cd->same_child && cd->same_child->nbestcolors) {
      if (!cd->nbestcolors || cd->same_child->nbestcolors < cd->nbestcolors) {
         int i;
         if (cd->nbestcolors) {
            COLORfree_sets(&(cd->bestcolors),&(cd->nbestcolors));
         }
         cd->nbestcolors = cd->same_child->nbestcolors;
         cd->same_child->nbestcolors = 0;
         cd->bestcolors = cd->same_child->bestcolors;
         cd->same_child->bestcolors = (COLORset*) NULL;
         for (i = 0; i < cd->nbestcolors; ++i) {
            int j; 
            int add_v2 = 0;
            for (j = 0; j < cd->bestcolors[i].count ; ++j) {
               if (cd->bestcolors[i].members[j] == cd->same_child->v1) {
                  add_v2 = 1;
               }
               if (cd->bestcolors[i].members[j] >= cd->same_child->v2) {
                  (cd->bestcolors[i].members[j])++;
               }
            }
            if (add_v2) {
               (cd->bestcolors[i].count)++;
               cd->bestcolors[i].members =
                  (int*) realloc(cd->bestcolors[i].members,
                                 cd->bestcolors[i].count * sizeof(int));
               COLORcheck_NULL(cd->bestcolors[i].members,
                               "Failed to realloc cd->bestcolors[i].members");
               cd->bestcolors[i].members[cd->bestcolors[i].count-1] = 
                  cd->same_child->v2;
            }
         }
         
         
      }
   }   
 CLEANUP:
   return rval;
}


static int collect_diff_child(colordata* cd)
{
   int rval = 0;

   if (cd->same_child && cd->diff_child->nbestcolors) {
      if (!cd->nbestcolors || cd->diff_child->nbestcolors < cd->nbestcolors) {
         if (cd->nbestcolors) {
            COLORfree_sets(&(cd->bestcolors),&(cd->nbestcolors));
         }
         cd->nbestcolors = cd->diff_child->nbestcolors;
         cd->diff_child->nbestcolors = 0;
         cd->bestcolors = cd->diff_child->bestcolors;
         cd->diff_child->bestcolors = (COLORset*) NULL;
      }
   }
   
   return rval;
}



static int compute_coloring(colordata* cd, int parent_lb, int parent_ub, int depth)
{
   int rval = 0;
   int iterations       = 0;
   int break_while_loop = 1;

   double last_snapshot_time;

   double cur_time;
   double lb_rtime;

   double child_lb = parent_lb; /*child_lb will be increased in this routine*/
   
   printf("Starting compute_coloring with lb %d and ub %d at depth %d.\n",
          parent_lb,parent_ub,depth);

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
      COLORlp_write(cd->lp,"debug.lp");
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

      if ((double) parent_lb < cd->dbl_lower_bound) {
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
             ((double) parent_lb < cd->dbl_lower_bound) && 
             !break_while_loop);

   if (iterations < cd->maxiterations) {
      char   mps_fname[256];
      printf ("Found bound of %lld (%20.16g), greedy coloring %d (iterations = %d).\n", 
              (long long) cd->lower_bound,cd->dbl_lower_bound, cd->ccount,iterations);
      lb_rtime = COLORcpu_time() - lb_rtime;
      printf("Computing initial lower bound took %f seconds.\n",lb_rtime);
      fflush(stdout);
      if (write_mwis) {
         sprintf(mps_fname,"%s.mwis.lp",cd->pname);
         COLORstable_write_mps(mps_fname,cd->ncount,cd->ecount,cd->elist,
                               cd->mwis_pi,cd->mwis_pi_scalef);

         sprintf(mps_fname,"%s.mwis.dimacs",cd->pname);
         rval = COLORstable_write_dimacs(mps_fname,cd->ncount,cd->ecount,cd->elist,
                                         cd->mwis_pi,cd->mwis_pi_scalef);

         sprintf(mps_fname,"%s.mwclq.dimacs",cd->pname);
         rval = COLORstable_write_dimacs_clique(mps_fname,
                                                cd->ncount,cd->ecount,cd->elist,
                                                cd->mwis_pi,cd->mwis_pi_scalef);

         sprintf(mps_fname,"%s.graphviz.dot",cd->pname);
         COLORplot_graphviz(mps_fname,cd->ncount,cd->ecount,cd->elist,(int*) NULL); 
      }
      if (cd->lower_bound < parent_ub) {
         /* Dive into child trees */

         delete_old_colorclasses(cd);
         
         if (cclasses_outfile) {
            int add_timestamp = 0;
            write_snapshot(cd,add_timestamp);
         }
         
         /*       { */
         /*          double colheur_rtime = COLORcpu_time(); */
         /*          heur_colors_with_stable_sets(cd); */
         /*          colheur_rtime = COLORcpu_time() - colheur_rtime; */
         /*          printf("Computing final upper bound took %f seconds.\n",colheur_rtime); */
         /*       } */
         
         /* Create branches or fill cd->bestcolors iff current LP
            solution is already integral
         */
         create_branches(cd);
         
         assert(cd->same_child || (cd->lower_bound == cd->nbestcolors));
                
         if (!cd->nbestcolors || cd->lower_bound < cd->nbestcolors) {
            if (cd->same_child) {
               printf("Entering same child\n");
               compute_coloring(cd->same_child, cd->lower_bound, parent_ub, depth+1);
               
               collect_same_child(cd);
               child_lb = cd->same_child->lower_bound;
               printf("Leaving same child\n");
            /*             free_colordata(cd->same_child);cd->same_child = (colordata*) NULL; */
            }
            
            if ( (cd->lower_bound < cd->nbestcolors) && cd->diff_child) {
               printf("Entering diff child\n");
               assert(cd->nbestcolors);
               compute_coloring(cd->diff_child,cd->lower_bound,cd->nbestcolors, depth+1);
            
               collect_diff_child(cd);
               
            if (cd->diff_child->lower_bound < child_lb) {
               child_lb = cd->diff_child->lower_bound;
            }
            printf("Leaving diff child\n");
            /*             free_colordata(cd->diff_child);cd->diff_child = (colordata*) NULL; */
            }
            assert(cd->lower_bound <= child_lb);
            if(cd->lower_bound <  child_lb) {
               cd->lower_bound = child_lb;
            }
         }
      }
   } else {
      printf ("Lower bound could not be found in %d iterations!\n", iterations);
   }
 CLEANUP:
   return rval;
}

int main (int ac, char **av)
{
    int rval = 0;
    double start_time;
    double tot_rtime;
   
    colordata  _colordata;
    colordata* cd = &_colordata;

    init_colordata(cd);
    
    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    get_problem_name(cd->pname,edgefile);


    if (debug) printf ("Debugging turned on\n");
    if (outfile) printf ("Output File: %s\n", outfile);
    fflush (stdout);


    rval = COLORread_dimacs (edgefile, &(cd->ncount), &(cd->ecount),
                             &(cd->elist), (int **) NULL);
    COLORcheck_rval (rval, "COLORread_diamcs failed");

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
    }

    compute_coloring(cd,1,cd->ccount,0);


    printf ("Opt Colors: %d\n", cd->nbestcolors); fflush (stdout);
    print_colors(cd->bestcolors,cd->nbestcolors);

    
    tot_rtime = COLORcpu_time() - start_time;
    printf("Computing coloring took %f seconds.\n",tot_rtime);
      

    
CLEANUP:
    free_colordata(cd);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int rval = 0;

    while ((c = getopt (ac, av, "dmo:r:w:")) != EOF) {
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
    fprintf (stderr, "   -r f  write initial stable sets from file f\n");
    fprintf (stderr, "   -w f  write stable sets from file f\n");
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

