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

#include "bbsafe.h"

#include "color_private.h"




static int debug = 0;
static int write_mwis = 0;
static int ncolordata = 0;

static char* cclasses_outfile;
static const int branch_with_same_sequence = 0;


typedef struct branching_joblist branching_joblist;
struct branching_joblist {
   colordata* cd;
   branching_joblist* next;
};


static int compute_lower_bound(colordata* cd);
static int grab_integral_solution(colordata* cd, double* x, double tolerance);

void COLORset_dbg_lvl(int dbglvl)
{
   debug = dbglvl;
}

void COLORset_cclasses_outfile(char* outfile)
{
   cclasses_outfile = outfile;
}

void COLORset_write_mwis(int write_mwis_par)
{
   write_mwis = write_mwis_par;
}


void init_colordata(colordata* cd)
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
   cd->upper_bound = COLORNWT_MAX;
   cd->lower_bound = 1;
   cd->dbl_safe_lower_bound = 1.0;
   cd->dbl_est_lower_bound = 1.0;
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

   cd->same_children = (colordata*)NULL;
   cd->nsame         = 0;
   cd->diff_children = (colordata*)NULL;
   cd->ndiff         = 0;
}

void free_lbcolordata(colordata* cd)
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
   cd->ccount = 0;

   COLORstable_freeenv(&(cd->mwis_env));

}

static 
void free_temporary_data(colordata* cd)
{
   int i;
   for (i = 0; i < cd->nsame;++ i) {
      free_colordata(&(cd->same_children[i]));
   }
   if (cd->same_children) free(cd->same_children);

   for (i = 0; i < cd->ndiff;++ i) {
      free_colordata(&(cd->diff_children[i]));
   }
   if (cd->diff_children) free(cd->diff_children);

   free_lbcolordata(cd);
   if (cd->elist) free (cd->elist);
   if (cd->orig_node_ids) free(cd->orig_node_ids);

}
void free_colordata(colordata* cd)
{
   free_temporary_data(cd);
   COLORfree_sets(&(cd->bestcolors),&(cd->nbestcolors));
}

COLOR_MAYBE_UNUSED
static int write_colordata(colordata* cd, FILE* file) {
   int rval = 0;
   int i;

   fprintf(file,"node %d\n",  cd->id);
   fprintf(file,"name %s\n",  cd->pname);
   fprintf(file,"depth %d\n", cd->depth);
   fprintf(file,"status %d\n", cd->status);
   fprintf(file,"ncount %d\n", cd->ncount);
   fprintf(file,"ecount %d\n", cd->ecount);

   fprintf(file,"edges");
   for (i = 0; i < cd->ecount; ++ i) {
      fprintf(file," %d %d", cd->elist[2*i],cd->elist[2*i+1]);
   }
   fprintf(file,"\n");

   fprintf(file,"orig_node_ids");
   for (i = 0; i < cd->ncount; ++ i) {
      fprintf(file,"orig_node_ids %d ",cd->orig_node_ids[i]);
   }
   fprintf(file,"\n");

   fprintf(file,"lower_bound %d\n",cd->lower_bound);
   fprintf(file,"upper_bound %d\n",cd->upper_bound);
   fprintf(file,"dbl_safe_lower_bound %f\n",cd->dbl_safe_lower_bound);

   fprintf(file,"ccount %d\n",cd->ccount);
   for (i = 0; i < cd->ccount; ++ i) {
      int j;
      fprintf(file,"cclass %d",cd->cclasses[i].count);
      for (j = 0; j < cd->cclasses[i].count; ++j) {
         fprintf(file," %d", cd->cclasses[i].members[j]);
      }
      fprintf(file,"\n");
   }


   fprintf(file,"nbestcolors %d\n", cd->nbestcolors);
   for (i = 0; i < cd->nbestcolors; ++ i) {
      int j;
      fprintf(file,"cclass %d",cd->bestcolors[i].count);
      for (j = 0; j < cd->bestcolors[i].count; ++j) {
         fprintf(file," %d", cd->bestcolors[i].members[j]);
      }
      fprintf(file,"\n");
   }

   fprintf(file,"nsame %d\n",cd->nsame);
   fprintf(file,"same_children");
   for (i = 0; i < cd->nsame; ++ i) {
      fprintf(file," %d",cd->same_children[i].id);
   }

   fprintf(file,"ndiff %d\n",cd->ndiff);
   fprintf(file,"diff_children");
   for (i = 0; i < cd->ndiff; ++ i) {
      fprintf(file," %d",cd->diff_children[i].id);
   }

   for (i = 0; i < cd->nsame; ++ i) {
      rval = write_colordata(cd->same_children + i,file);
      COLORcheck_rval(rval,"Failed in write_colordata");
   }


   for (i = 0; i < cd->ndiff; ++ i) {
      rval = write_colordata(cd->diff_children + i,file);
      COLORcheck_rval(rval,"Failed in write_colordata");
   }

 CLEANUP:
   return rval;
}

int send_colordata(COLOR_SFILE *s, colordata* cd, int include_best) {
   int rval = 0;
   int i;


   rval = COLORsafe_swrite_int (s,cd->id);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int (s,cd->id)");

   rval = COLORsafe_swrite_string (s,cd->pname);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_string(s,cd->pname)");

   rval = COLORsafe_swrite_int (s,cd->depth);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->depth)");

   rval = COLORsafe_swrite_uint(s,cd->status);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_uint(s,cd->status)");

   rval = COLORsafe_swrite_int(s,cd->ncount);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->ncount)");

   rval = COLORsafe_swrite_int(s,cd->ecount);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->ecount)");

   for (i = 0; i < cd->ecount; ++ i) {
      rval = COLORsafe_swrite_int(s,cd->elist[2*i]);
      COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->elist[2*i])");

      rval = COLORsafe_swrite_int(s,cd->elist[2*i + 1]);
      COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->elist[2*i + 1])");
   }

   for (i = 0; i < cd->ncount; ++ i) {
      rval = COLORsafe_swrite_int(s,cd->orig_node_ids[i]);
      COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->orig_node_ids[i])");
   }

   rval = COLORsafe_swrite_int(s,cd->lower_bound);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->lower_bound)");

   rval = COLORsafe_swrite_int(s,cd->upper_bound);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->upper_bound)");

   rval = COLORsafe_swrite_double(s,cd->dbl_safe_lower_bound);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_double(s,cd->dbl_safe_lower_bound)");

   rval = COLORsafe_swrite_double(s,cd->dbl_est_lower_bound);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_double(s,cd->dbl_est_lower_bound)");

   rval = COLORsafe_swrite_int(s,cd->ccount);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->ccount)");
   for (i = 0; i < cd->ccount; ++ i) {
      int j;
      rval = COLORsafe_swrite_int(s,cd->cclasses[i].count);
      COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->cclasses[i].count)");
      for (j = 0; j < cd->cclasses[i].count; ++j) {
         rval = COLORsafe_swrite_int(s,cd->cclasses[i].members[j]);
         COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->cclasses[i].members[j])");
      }
   }


   if (include_best) {
      rval = COLORsafe_swrite_int(s,cd->nbestcolors);
      COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->nbestcolors)");

      for (i = 0; i < cd->nbestcolors; ++ i) {
         int j;
         rval = COLORsafe_swrite_int(s,cd->bestcolors[i].count);
         COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->bestcolors[i].count)");
         for (j = 0; j < cd->bestcolors[i].count; ++j) {
            rval = COLORsafe_swrite_int(s,cd->bestcolors[i].members[j]);
            COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->bestcolors[i].members[j])");
         }
      }
   }


   rval = COLORsafe_swrite_int(s,cd->v1);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->v1)");

   rval = COLORsafe_swrite_int(s,cd->v2);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->v2)");



   rval = COLORsafe_swrite_int(s,cd->nsame);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->nsame)");

   rval = COLORsafe_swrite_int(s,cd->ndiff);
   COLORcheck_rval(rval,"Failed in COLORsafe_swrite_int(s,cd->ndiff)");

   for (i = 0; i < cd->nsame; ++ i) {
      rval = send_colordata(s,cd->same_children + i,include_best);
      COLORcheck_rval(rval,"Failed in sent_colordata");
   }

   for (i = 0; i < cd->ndiff; ++ i) {
      rval = send_colordata(s,cd->diff_children + i,include_best);
      COLORcheck_rval(rval,"Failed in sent_colordata");
   }

 CLEANUP:
   return rval;
}

static
int receive_child_data(COLOR_SFILE *s,
                       colordata** children,int nchildren,
                       colordata* parent,
                       int adopt_id, int include_best)
{
   int rval = 0;
   int i;
   if (nchildren) {
      *children = COLOR_SAFE_MALLOC(nchildren,colordata);
      
      for (i = 0; i < nchildren; ++ i) {
         init_colordata(*children + i);
         rval = receive_colordata(s,*children + i,adopt_id,include_best);
         if (!adopt_id) {
            (*children)[i].parent = parent;
         }
         COLORcheck_rval(rval,"Failed in sent_colordata");
      }
   } else {
      *children = (colordata*) NULL;
   }
 CLEANUP:
   return rval;
}

int receive_colordata(COLOR_SFILE *s, colordata* cd, 
                      int adopt_id,
                      int include_best) {
   int rval = 0;
   int i;
   int new_id;

   rval = COLORsafe_sread_int (s,&new_id);
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_int (s,cd->id)");
   if (adopt_id && new_id != cd->id) {
      cd->id= new_id;
   }

   rval = COLORsafe_sread_string(s,cd->pname,MAX_PNAME_LEN - 1);
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_string(s,cd->pname)");

   rval = COLORsafe_sread_int(s,&(cd->depth));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->depth)");

   rval = COLORsafe_sread_uint(s, &(cd->status));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_uint(s,cd->status)");

   rval = COLORsafe_sread_int(s,&(cd->ncount));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->ncount)");

   rval = COLORsafe_sread_int(s,&(cd->ecount));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->ecount)");

   cd->elist = (int*) COLOR_SAFE_MALLOC(2 * cd->ecount, int);
   COLORcheck_NULL(cd->elist, "Failed to allocate cd->elist");

   for (i = 0; i < cd->ecount; ++ i) {
      rval = COLORsafe_sread_int(s,cd->elist + 2*i);
      COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->elist + 2*i)");

      rval = COLORsafe_sread_int(s,cd->elist + 2*i + 1);
      COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->elist + 2*i + 1)");
   }

   cd->orig_node_ids = COLOR_SAFE_MALLOC(cd->ncount,int);
   COLORcheck_NULL(cd->orig_node_ids,"Failed to allocate cd->orig_node_ids");

   for (i = 0; i < cd->ncount; ++ i) {
      rval = COLORsafe_sread_int(s,cd->orig_node_ids + i);
      COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->orig_node_ids + i)");
   }

   rval = COLORsafe_sread_int(s,&(cd->lower_bound));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->lower_bound)");

   rval = COLORsafe_sread_int(s,&(cd->upper_bound));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->upper_bound)");

   rval = COLORsafe_sread_double(s,&(cd->dbl_safe_lower_bound));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_double(s,cd->dbl_safe_lower_bound)");

   rval = COLORsafe_sread_double(s,&(cd->dbl_est_lower_bound));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_double(s,cd->dbl_est_lower_bound)");

   rval = COLORsafe_sread_int(s,&(cd->ccount));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,&(cd->ccount))");
   cd->gallocated = cd->ccount;
   if (cd->cclasses) {
      free (cd->cclasses);
      cd->cclasses = (COLORset*) NULL;
   }
   if (cd->ccount) {
      assert(!cd->cclasses);
      cd->cclasses = COLOR_SAFE_MALLOC(cd->ccount, COLORset);
      COLORcheck_NULL(cd->cclasses,"Failed to allocate cd->cclasses");
      for (i = 0; i < cd->ccount; ++ i) {
         int j;
         rval = COLORsafe_sread_int(s,&(cd->cclasses[i].count));
         COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->cclasses[i].count)");

         cd->cclasses[i].members = COLOR_SAFE_MALLOC(cd->cclasses[i].count,int);
         COLORcheck_NULL(cd->cclasses[i].members,"Failed to allocate cd->cclasses[i].members");
         
         for (j = 0; j < cd->cclasses[i].count; ++j) {
            rval = COLORsafe_sread_int(s,&(cd->cclasses[i].members[j]));
            COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->cclasses[i].members[j])");
         }
      }
   }

   if (include_best) {
      int nbest;
      rval = COLORsafe_sread_int(s,&(nbest));
      COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->nbestcolors)");

      if (nbest) {
         cd->nbestcolors = nbest;

         COLOR_IFFREE(cd->bestcolors, COLORset);
         cd->bestcolors = COLOR_SAFE_MALLOC(cd->nbestcolors,COLORset);
         COLORcheck_NULL(cd->bestcolors,"Failed to allocate cd->bestcolors");
     
         
         for (i = 0; i < cd->nbestcolors; ++ i) {
            int j;
            rval = COLORsafe_sread_int(s,&(cd->bestcolors[i].count));
            COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->bestcolors[i].count)");
         
            cd->bestcolors[i].members = COLOR_SAFE_MALLOC(cd->bestcolors[i].count,int);
            COLORcheck_NULL(cd->bestcolors[i].members,"Failed to allocate cd->bestcolors[i].members");
         
            for (j = 0; j < cd->bestcolors[i].count; ++j) {
               rval = COLORsafe_sread_int(s,&(cd->bestcolors[i].members[j]));
               COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->bestcolors[i].members[j])");
            }
         }
      }
   }

   rval = COLORsafe_sread_int(s,&(cd->v1));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->v1)");

   rval = COLORsafe_sread_int(s,&(cd->v2));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->v2)");
   
   
   rval = COLORsafe_sread_int(s,&(cd->nsame));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->nsame)");

   rval = COLORsafe_sread_int(s,&(cd->ndiff));
   COLORcheck_rval(rval,"Failed in COLORsafe_sread_int(s,cd->ndiff)");


   rval = receive_child_data(s,&(cd->same_children), cd->nsame,
                             cd,
                             adopt_id,include_best);
   COLORcheck_rval(rval,"Failed in receive_child_data");

   rval = receive_child_data(s,&(cd->diff_children), cd->ndiff,
                             cd,
                             adopt_id,include_best);
   COLORcheck_rval(rval,"Failed in receive_child_data");

/*    if (cd->nsame) { */
/*       cd->same_children = COLOR_SAFE_MALLOC(cd->nsame,colordata); */

/*       for (i = 0; i < cd->nsame; ++ i) { */
/*          init_colordata(cd->same_children + i); */
/*          rval = receive_colordata(s,cd->same_children + i,adopt_id,include_best); */
/*          if (!adopt_id) { */
/*             cd->same_children[i].parent = cd; */
/*          } */
/*          COLORcheck_rval(rval,"Failed in sent_colordata"); */
/*       } */
/*    } else { */
/*       cd->same_children = (colordata*) NULL; */
/*    } */

/*    if (cd->ndiff) { */
/*       cd->diff_children = COLOR_SAFE_MALLOC(cd->ndiff,colordata); */

/*       for (i = 0; i < cd->ndiff; ++ i) { */
/*          init_colordata(cd->diff_children + i); */
/*          rval = receive_colordata(s,cd->diff_children + i,adopt_id,include_best); */
/*          COLORcheck_rval(rval,"Failed in sent_colordata"); */
/*          if (!adopt_id) { */
/*             cd->same_children[i].parent = cd; */
/*          } */
/*       } */
/*    } else { */
/*       cd->diff_children = (colordata*) NULL; */
/*    } */

 CLEANUP:
   return rval;
}


int COLORdbg_lvl() {
   return debug;
}

int print_colors(COLORset* cclasses, int ccount)
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

   cd->dbl_safe_lower_bound = COLORsafe_lower_dbl(cd->lower_scaled_bound,cd->mwis_pi_scalef);
   if (COLORdbg_lvl() > 0) {
      printf("Current primal LP objective: %f (%lld / %lld).\n",
             cd->dbl_safe_lower_bound,
             (long long) cd->lower_scaled_bound,
             (long long) cd->mwis_pi_scalef );
   }
   cd->lower_bound = ceil(cd->dbl_safe_lower_bound);

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

   if (numdel) {
   
      if (COLORdbg_lvl() > 0) {
         printf("Deleted %d out of %d columns with age > %d. Rebuilding LP from scratch.\n",
                numdel, numdel + cd->ccount, cd->retirementage);
      }
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

int build_lp(colordata* cd)
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

static int mark_neighborhood(int* neighbor_marker,
                             int ncount, int ecount, int elist[],
                             int v)

{
   int rval = 0;
   int i;
   COLORadjgraph G;

   COLORadjgraph_init(&G);
   rval = COLORadjgraph_build(&G,ncount,ecount, elist);
   COLORcheck_rval(rval,"COLORadjgraph_build");

   for(i = 0; i < ncount; ++i) { neighbor_marker[i] = 0;}

   for(i = 0; i < G.nodelist[v].degree; ++i) {
      int v_i = G.nodelist[v].adj[i];
      neighbor_marker[v_i] = 1;
   }
   COLORadjgraph_free(&G);

 CLEANUP:
   return rval;
}

static int create_contracted_graph(colordata* cd,
                                   int ncount, int ecount, int elist[],
                                   int v1, int v2)
{
   int rval = 0;
   int i;
   COLORadjgraph G;




   /* Create contracted graph */
   cd->ncount = ncount - 1;
   cd->ecount = ecount;
   cd->elist  = (int*) COLOR_SAFE_MALLOC(2 * cd->ecount,int);
   COLORcheck_NULL(cd->elist,"Failed to allocate cd->elist");
   memcpy(cd->elist,elist, 2 * ecount * sizeof(int));
   for(i = 0; i < cd->ecount; ++i) {
      cd->elist[2*i] = new_eindex(cd->elist[2*i],v1,v2);
      cd->elist[2*i+1] = new_eindex(cd->elist[2*i+1],v1,v2);
   }

   COLORadjgraph_init(&G);
   rval = COLORadjgraph_build(&G,cd->ncount,cd->ecount, cd->elist);
   COLORcheck_rval(rval,"COLORadjgraph_build");

   rval =  COLORadjgraph_simplify(&G);
   COLORcheck_rval(rval,"COLORadjgraph_simplify");

   rval = COLORadjgraph_extract_edgelist(&(cd->ecount), &(cd->elist),&G);

   cd->orig_node_ids = (int*) COLOR_SAFE_MALLOC(cd->ncount,int);
   COLORcheck_NULL(cd->orig_node_ids,"Failed to allocate cd->orig_node_ids");
   for (i = 0; i < cd->ncount; ++i) {
      int j = (i < v2 ) ?  i : i + 1;
      cd->orig_node_ids[i] = cd->parent->orig_node_ids[j];
   }

 CLEANUP:
   COLORadjgraph_free(&G);

   return rval;
}

static int prune_duplicated_sets (colordata* cd)
{
   int rval = 0;
   int i,j;
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
 CLEANUP:
   return rval;
}

static int transfer_same_cclasses(colordata* cd,
                                  const int* v2_neighbor_marker,
                                  const COLORset* parent_cclasses,
                                  int   parent_ccount,
                                  int   v1,
                                  int   v2)
{
   int rval = 0;
   int i;
   int v1_covered = 0;
   /* Transfer independent sets: */
   cd->gallocated = cd->ccount   =  parent_ccount + 1;
   cd->cclasses = (COLORset*) COLOR_SAFE_MALLOC(cd->gallocated,COLORset);
   for (i = 0; i < parent_ccount; ++i) {
      int j;
      int add_v1 = 1;
      cd->cclasses[i].members = (int*) COLOR_SAFE_MALLOC(parent_cclasses[i].count,int);
      cd->cclasses[i].count = 0;
      for (j = 0; j < parent_cclasses[i].count; ++j) {
         if (v2_neighbor_marker[parent_cclasses[i].members[j]] == 1) {
            add_v1 = 0;
            j = parent_cclasses[i].count;/*break*/
         }
      }
      for (j = 0; j < parent_cclasses[i].count; ++j) {
         if (parent_cclasses[i].members[j] == v1) {
            if (add_v1) {v1_covered = 1;}
            else {continue;}
         }
         if (parent_cclasses[i].members[j] < v2) {
            cd->cclasses[i].members[(cd->cclasses[i].count)++] =
               parent_cclasses[i].members[j];
         }
         else if (parent_cclasses[i].members[j] >  v2) {
            cd->cclasses[i].members[(cd->cclasses[i].count)++] =
               parent_cclasses[i].members[j] - 1;
         }
         /* else 'parent_cclasses[i].members[j] == v2' and we skip it*/
      }
      if(COLORdbg_lvl() > 1) {
         printf("PARENT SET ");
         for (j = 0; j < parent_cclasses[i].count; ++j) {
            printf(" %d",parent_cclasses[i].members[j]);
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
      cd->cclasses[parent_ccount].count  = 1;
      cd->cclasses[parent_ccount].members = (int*) COLOR_SAFE_MALLOC(1,int);
      cd->cclasses[parent_ccount].members[0] = v1;
   } else {
      cd->cclasses[parent_ccount].count  = 0;
      cd->cclasses[parent_ccount].members = (int*) 0;
      cd->ccount--;
   }
   rval = prune_duplicated_sets(cd);
   COLORcheck_rval(rval,"Failed in prune_duplicated_sets");
   /* END Transfer independent sets: */

 CLEANUP:
   return rval;
}

static int create_same (colordata* parent_cd, int v1, int v2)
{
   int rval = 0;
   colordata*    cd = (colordata*) NULL;
   int* v2_neighbor_marker = (int*) COLOR_SAFE_MALLOC(parent_cd->ncount,int);
   COLORcheck_NULL(v2_neighbor_marker,"Failed to allocate v2_neighbor_marker");

   cd = (colordata*) COLOR_SAFE_MALLOC(1,colordata);
   COLORcheck_NULL(cd,"Failed to allocate cd");

   init_colordata(cd);
   cd->depth = parent_cd->depth + 1;
   parent_cd->nsame         = 1;
   parent_cd->same_children = cd;

   cd->v1 = v1;
   cd->v2 = v2;

   cd->upper_bound = parent_cd->upper_bound;
   cd->lower_bound = parent_cd->lower_bound;
   cd->dbl_safe_lower_bound = parent_cd->dbl_safe_lower_bound;

   cd->parent = parent_cd;
   cd->debugcolors = parent_cd->debugcolors;
   cd->ndebugcolors = parent_cd->ndebugcolors;

   rval = mark_neighborhood(v2_neighbor_marker,
                            parent_cd->ncount,parent_cd->ecount,parent_cd->elist,
                            v2);
   COLORcheck_rval(rval,"Failed in mark_neighborhood");

   rval = create_contracted_graph (cd,
                                   parent_cd->ncount,
                                   parent_cd->ecount,parent_cd->elist,
                                   v1,v2);
   COLORcheck_rval(rval,"Failed in create_contracted_graph_and_intersect");


   /* END create contracted graph */

   if (COLORdbg_lvl() > 1) {
      printf("create_same created following graph:\n");
      COLORgraph_print(cd->ecount,cd->elist);
   }

   rval = transfer_same_cclasses(cd,
                                 v2_neighbor_marker,
                                 parent_cd->cclasses,
                                  parent_cd->ccount,
                                 v1,v2);
   COLORcheck_rval(rval,"Failed in transfer_same_cclasses");
 CLEANUP:
   if (rval) {
      if (cd) {
         free_colordata(cd);
         free(cd);
      }
      parent_cd->same_children = (colordata*) NULL;
   }
   if (v2_neighbor_marker) free(v2_neighbor_marker);
   return rval;
}


static int create_differ(colordata* parent_cd, int v1, int v2)
{
   int rval = 0;
   int i;
   COLORadjgraph G;
   int           v2_covered = 0;
   colordata*    cd = (colordata*) COLOR_SAFE_MALLOC(1,colordata);
   COLORcheck_NULL(cd,"Failed to allocate cd");


   init_colordata(cd);

   cd->depth = parent_cd->depth + 1;
   parent_cd->ndiff         = 1;
   parent_cd->diff_children = cd;

   cd->v1 = v1;
   cd->v2 = v2;

   cd->upper_bound = parent_cd->upper_bound;
   cd->lower_bound = parent_cd->lower_bound;
   cd->dbl_safe_lower_bound = parent_cd->dbl_safe_lower_bound;


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


   rval = prune_duplicated_sets(cd);
   COLORcheck_rval(rval,"Failed in prune_duplicated_sets");


 CLEANUP:
   if (rval) {
      if (cd) {
         free_colordata(cd);
         free(cd);
      }
      parent_cd->diff_children = (colordata*) NULL;
   }
   COLORadjgraph_free(&G);

   return rval;
}


static int apply_diffpair_to_cclasses(COLORset** cclasses_ptr,
                                      int*       ccount_ptr,
                                      int        v1,
                                      int        v2)
{
   int rval = 0;
   int v1_covered = 0;
   int v2_covered = 0;
   int i;
   COLORset* cclasses = *cclasses_ptr;
   int       ccount   = *ccount_ptr;

   for (i = 0; i < ccount; ++i) {
      int v1_in_class = 0;
      int j;
      for (j = 0; j < cclasses[i].count; ++j) {
         int current_elm = cclasses[i].members[j];
         if (current_elm == v1) {
            if (v1_covered && !v2_covered) {
               memmove( cclasses[i].members + j,
                        cclasses[i].members + j + 1,
                        (cclasses[i].count - 1 - j)* sizeof(int));
               cclasses[i].count--;
            } else {
               v1_in_class = 1;
               v1_covered = 1;
            }
         } else if (current_elm == v2) {
            if (v1_in_class) {
               memmove( cclasses[i].members + j,
                        cclasses[i].members + j + 1,
                        (cclasses[i].count - 1 - j)* sizeof(int));
               cclasses[i].count--;
            } else {
               v2_covered = 1;
            }
         }
      }
   }
   if (!v2_covered) {
      cclasses = (COLORset*) realloc(cclasses, ccount + 1);
      COLORcheck_NULL(cclasses,"Failed to realloc cclasses");
      cclasses[ccount].members = (int*) COLOR_SAFE_MALLOC(1,int);
      cclasses[ccount++].members[0] = v2;
      cclasses[ccount++].count      = 1;
   }

   *cclasses_ptr = cclasses;
   *ccount_ptr   = ccount;
 CLEANUP:
   return rval;
}

/* Alternatively to the same and diff branch, we
   can create a same branch with v2 and each neighbor of v1>
*/
static int create_same_seq (colordata* parent_cd, int v1, int v2)
{
   int rval = 0;
   COLORadjgraph G;
   COLORadjnode* v1_node;
   int           i;
   colordata*    cd = (colordata*) NULL;
   int           nsame;
   int           ecount  = parent_cd->ecount;
   int*          tmp_elist = (int*) NULL;
   COLORset*     tmp_cclasses;
   int           tmp_ccount;
   int*          v2_neighbor_marker = (int*) NULL;
   int*          parent_v2_neighbor_marker =
      (int*) COLOR_SAFE_MALLOC(parent_cd->ncount,int);
   COLORcheck_NULL(parent_v2_neighbor_marker,"Failed to allocate parent_v2_neighbor_marker");

   v2_neighbor_marker =
      (int*) COLOR_SAFE_MALLOC(parent_cd->ncount,int);
   COLORcheck_NULL(v2_neighbor_marker,"Failed to allocate v2_neighbor_marker");

   rval = COLORcopy_sets(&tmp_cclasses,&tmp_ccount,
                         parent_cd->cclasses,parent_cd->ccount);
   COLORcheck_rval(rval,"Failed in COLORcopy_sets");

   COLORadjgraph_init(&G);
   rval = COLORadjgraph_build(&G,parent_cd->ncount,parent_cd->ecount, parent_cd->elist);
   COLORcheck_rval(rval,"COLORadjgraph_build");
   v1_node = &(G.nodelist[v1]);
   parent_cd->same_children = cd;


   tmp_elist = (int*) COLOR_SAFE_MALLOC(2*(parent_cd->ecount+v1_node->degree-1),
                                        int);
   COLORcheck_NULL(tmp_elist,"Failed to allocate cd");
   memcpy (tmp_elist, parent_cd->elist, 2 * parent_cd->ecount * sizeof(int));

   rval = create_same(parent_cd,v1,v2);
   COLORcheck_rval(rval,"Failed in create_same");
   nsame = 1;

   cd = (colordata*) COLOR_SAFE_MALLOC(v1_node->degree + 1,
                                       colordata);
   COLORcheck_NULL(cd, "Failed to allocate cd");
   for (i = 0; i < v1_node->degree + 1; ++i) {init_colordata(cd + i);}


   memcpy(cd,parent_cd->same_children,sizeof(colordata));
   free(parent_cd->same_children);
   parent_cd->same_children = cd;


   mark_neighborhood(parent_v2_neighbor_marker,
                     parent_cd->ncount,parent_cd->ecount,parent_cd->elist,v2);
   for (i = 0; i < v1_node->degree; ++i) {
      int n1_i = v1_node->adj[i];

      /** if n1_i is neighbor of v2, they cannot be colored equally.*/
      if (parent_v2_neighbor_marker[n1_i] == 1) {
         continue;
      }

      cd = &(parent_cd->same_children[nsame++]);

      cd->depth = parent_cd->depth + 1;
      cd->v1 = n1_i < v2 ? n1_i : v2;
      cd->v2 = n1_i > v2 ? n1_i : v2;

      if (COLORdbg_lvl() > 1) {
         printf("Creating same for %d %d.\n",cd->v1,cd->v2);
      }

      cd->upper_bound = parent_cd->upper_bound;
      cd->lower_bound = parent_cd->lower_bound;
      cd->dbl_safe_lower_bound = parent_cd->dbl_safe_lower_bound;

      cd->parent = parent_cd;
      cd->debugcolors = parent_cd->debugcolors;
      cd->ndebugcolors = parent_cd->ndebugcolors;

      mark_neighborhood(v2_neighbor_marker,
                        parent_cd->ncount,ecount,tmp_elist,cd->v2);

      rval = create_contracted_graph (cd,
                                      parent_cd->ncount,ecount,tmp_elist,
                                      cd->v1,cd->v2);
      COLORcheck_rval(rval,"create_contracted_graph_and_intersect");

      rval = transfer_same_cclasses(cd,
                                    v2_neighbor_marker,
                                    tmp_cclasses,
                                    tmp_ccount,
                                    cd->v1,cd->v2);
      COLORcheck_rval(rval,"Failed in transfer_same_cclasses");

      if (i + 1 > v1_node->degree) {
         /** Insert edge v1,n1_i to create DIFFER(v1,n1_i).*/
         tmp_elist[2 * ecount]         = v1;
         tmp_elist[2 * ecount + 1]     = n1_i;
         ecount++;
         rval = apply_diffpair_to_cclasses(&tmp_cclasses,&tmp_ccount,n1_i,v2);
         COLORcheck_rval(rval,"Failed in apply_diffpair_to_cclasses");
      }

      rval = compute_lower_bound(cd);
      COLORcheck_rval(rval,"Failed in compute_lower_bound");
   }

   parent_cd->nsame         = nsame;
   assert(parent_cd->ndiff == 0);
   assert(parent_cd->nsame >  0);

 CLEANUP:
   COLORadjgraph_free(&G);
   if (v2_neighbor_marker) free(v2_neighbor_marker);
   if (parent_v2_neighbor_marker) free(parent_v2_neighbor_marker);
   if (tmp_elist)  free(tmp_elist);
   COLORfree_sets(&tmp_cclasses,&tmp_ccount);
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
   assert(cd->nbestcolors == (int) ceil(incumbent));

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

   int    max_non_improving_branches  = 3; /* cd->ncount / 100 + 1; */
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
      if (COLORdbg_lvl() > 1) {
         printf("Creating branches for v1 = %d, v2 = %d (node-pair weight %f)\n",v1,v2,
                nodepair_weights[(int)(min_nodepair-nodepair_refs)]);
      }
      /* Create DIFFER and SAME */
      rval = create_same (cd,v1,v2);
      COLORcheck_rval(rval, "Failed in create_same");

      rval = create_differ(cd,v1,v2);
      COLORcheck_rval(rval, "Failed in create_differ");

      cd->same_children->maxiterations = 5;
      cd->diff_children->maxiterations = 5;

      compute_lower_bound(cd->same_children);
      compute_lower_bound(cd->diff_children);

      dbl_child_lb = (cd->same_children->dbl_safe_lower_bound < cd->diff_children->dbl_safe_lower_bound) ?
         cd->same_children->dbl_safe_lower_bound : cd->diff_children->dbl_safe_lower_bound;

      free_colordata(cd->same_children);
      cd->nsame=0; free(cd->same_children); cd->same_children=(colordata*) NULL;
      free_colordata(cd->diff_children);
      cd->ndiff=0;free(cd->diff_children); cd->diff_children = (colordata*) NULL;

      if (dbl_child_lb > strongest_dbl_lb) {
         strongest_dbl_lb = dbl_child_lb;
         *strongest_v1     = v1;
         *strongest_v2     = v2;
         remaining_branches = max_non_improving_branches;
      }
      if (COLORdbg_lvl() > 1) {
         printf("Found child bound of %f for v1 = %d, v2 = %d, nodepair_weight = %f .\n",dbl_child_lb,
                v1, v2,
                nodepair_weights[(int)(min_nodepair-nodepair_refs)]);
      }
   }
   {
      int nodepair_ref = nodepair_ref_key(*strongest_v1,*strongest_v2);

      if (COLORdbg_lvl() > 0) {
         printf("Found strongest child bound of %f for v1 = %d, "
                "v2 = %d, nodepair_weight = %f .\n",
                strongest_dbl_lb, *strongest_v1,
                *strongest_v2, nodepair_weights[nodepair_ref]);
      }
   }
 CLEANUP:
   return rval;
}

int create_branches(colordata* cd)
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

         
         
   if (COLORdbg_lvl() > 1) {
      printf("Collected %d branching candidates.\n",COLORNWTheap_size(cand_heap));
   }

   rval = find_strongest_children(&strongest_v1,&strongest_v2,
                                  cd,cand_heap,nodepair_refs,
                                  nodepair_weights);
   COLORcheck_rval(rval, "Failed in find_strongest_children");


   /*    Create DIFFER and SAME for strongest children */
   if (branch_with_same_sequence) {
      rval = create_same_seq(cd,strongest_v1,strongest_v2);
      COLORcheck_rval(rval, "Failed in create_same_seq");
   } else { /* branch with one same and one diff branch*/
      rval = create_same (cd,strongest_v1,strongest_v2);
      COLORcheck_rval(rval, "Failed in create_same");
      rval = compute_lower_bound (cd->same_children);
      COLORcheck_rval(rval, "Failed in compute_lower_bound");

      rval = create_differ(cd,strongest_v1,strongest_v2);
      COLORcheck_rval(rval, "Failed in create_differ");

      rval = compute_lower_bound (cd->diff_children);
      COLORcheck_rval(rval, "Failed in compute_lower_bound");
   }

   if (cd->same_children && cd->ndebugcolors) {
      int same_opt_track = 0;

      if (cd->opt_track) {
         same_opt_track = are_in_same(cd->debugcolors,cd->ndebugcolors,
                                      cd->orig_node_ids[cd->same_children->v1],
                                      cd->orig_node_ids[cd->same_children->v2]);
         cd->same_children->opt_track = same_opt_track;
         cd->diff_children->opt_track = !same_opt_track;
      } else  {
         cd->same_children->opt_track = 0;
         cd->diff_children->opt_track = 0;
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

/** Transform the coloring of cd->same_children into a coloring of cd.
 */
static int collect_same_children(colordata* cd)
{
   int rval = 0;
   int c;

   for (c = 0; c < cd->nsame; ++c) {
      if ( cd->same_children[c].nbestcolors &&
           (!cd->nbestcolors || cd->same_children[c].nbestcolors < cd->upper_bound) ) {
         int i;
         if (cd->nbestcolors) {
            COLORfree_sets(&(cd->bestcolors),&(cd->nbestcolors));
         }
         cd->upper_bound = cd->nbestcolors = cd->same_children[c].nbestcolors;
         cd->same_children[c].nbestcolors = 0;
         cd->bestcolors = cd->same_children[c].bestcolors;
         cd->same_children[c].bestcolors = (COLORset*) NULL;
         for (i = 0; i < cd->nbestcolors; ++i) {
            int j;
            int add_v2 = 0;
            for (j = 0; j < cd->bestcolors[i].count ; ++j) {
               if (cd->bestcolors[i].members[j] == cd->same_children[c].v1) {
                  assert(add_v2 == 0); add_v2 = 1;
               }
               if (cd->bestcolors[i].members[j] >= cd->same_children[c].v2) {
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
                       (cd->same_children[c].v2 < cd->bestcolors[i].members[k - 1]) ) {
                  k--;
               }
               if (k < cd->bestcolors[i].count - 1) {
                  memmove(cd->bestcolors[i].members + k + 1,
                          cd->bestcolors[i].members + k,
                          (cd->bestcolors[i].count - 1 - k) * sizeof(int));
               }
               cd->bestcolors[i].members[k] = cd->same_children[c].v2;

            }
         }
/*          print_colors(cd->bestcolors,cd->nbestcolors); */
         fflush(stdout);
         rval = COLORcheck_coloring(cd->bestcolors,cd->nbestcolors,
                                    cd->ncount, cd->ecount, cd->elist);
         COLORcheck_rval(rval,"ERROR: An incorrect coloring was created.");
      }

   }
 CLEANUP:
   return rval;
}

/** Transform the coloring of cd->same_children into a coloring of cd.
 */
static int collect_diff_children(colordata* cd)
{
   int rval = 0;
   int c;

   for (c = 0; c < cd->ndiff; ++c) {
      if ( cd->diff_children[c].nbestcolors &&
           ( !cd->nbestcolors ||
             cd->diff_children[c].nbestcolors < cd->upper_bound) )
         {
            if (cd->nbestcolors) {
               COLORfree_sets(&(cd->bestcolors),&(cd->nbestcolors));
            }
            cd->upper_bound = cd->nbestcolors = cd->diff_children[c].nbestcolors;
            cd->diff_children[c].nbestcolors = 0;
            cd->bestcolors = cd->diff_children[c].bestcolors;
            cd->diff_children[c].bestcolors = (COLORset*) NULL;

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

   if (COLORdbg_lvl() > 1) {
      printf("Starting compute_lower_bound with lb %d (%f) and ub %d at depth %d (id = %d, opt_track = %d).\n",
             cd->lower_bound,cd->dbl_safe_lower_bound,cd->upper_bound,cd->depth,cd->id, cd->opt_track);
   }

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
         rval = delete_old_colorclasses(cd);
         COLORcheck_rval (rval, "delete_old_colorclasses failed");
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
/*       printf("COLORlp_optimize succeeded\n");fflush(stdout); */
      rval = grow_ages(cd);
      COLORcheck_rval (rval, "grow_ages failed");

      if (COLORdbg_lvl() > 1) {print_ages(cd);}

      rval = COLORlp_pi (cd->lp, cd->pi);
      COLORcheck_rval (rval, "COLORlp_pi failed");

      make_pi_feasible(cd);

      COLOR_double2COLORNWT(cd->mwis_pi,&(cd->mwis_pi_scalef),
                            cd->pi,cd->ncount);

      compute_objective(cd);
      cd->dbl_est_lower_bound = cd->dbl_safe_lower_bound;

      if (iterations < cd->maxiterations) {
         int set_i;

         /** Solve the MWIS problem. The node weigths might be scaled
             down further in this method!*/
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
             !break_while_loop);
   if (iterations < cd->maxiterations) {
      /** As the node weigths can have been scaled down further since
          the last computation of the objective, we need to recompute
          it here.*/
      compute_objective(cd);

      if (COLORdbg_lvl() > 1) {
         printf ("Found bound of %lld (%20.16g), upper_bound %d (id = %d, iterations = %d, opt_track = %d).\n",
                 (long long) cd->lower_bound,cd->dbl_safe_lower_bound,
                 cd->upper_bound,cd->id,iterations, cd->opt_track);
      }

      cd->retirementage = 1;
      rval = delete_old_colorclasses(cd);
      COLORcheck_rval (rval, "delete_old_colorclasses failed");

      lb_rtime = COLORcpu_time() - lb_rtime;
      cd->status = LP_bound_computed;
      if (COLORdbg_lvl()>=0) {
         printf("Computing initial lower bound took %f seconds.\n",lb_rtime);
      }
   } else {
      lb_rtime = COLORcpu_time() - lb_rtime;
      cd->status = LP_bound_estimated;
      if (COLORdbg_lvl()>=0) {
         printf("Computing estimated lower bound took %f seconds.\n",lb_rtime);
      }
   }
   fflush(stdout);

 CLEANUP:
   return rval;
}

static int trigger_lb_changes(colordata* child)
{
   int rval = 0;
   int i;
   int new_lower_bound = child->lower_bound;

   colordata* cd = (colordata*) child->parent;
   
   while (cd) {
      for (i = 0;  i < cd->nsame; ++i) {
         if (cd->same_children[i].lower_bound < new_lower_bound) {
            new_lower_bound = cd->same_children[i].lower_bound;
         }
      }
      for (i = 0;  i < cd->ndiff; ++i) {
         if (cd->diff_children[i].lower_bound < new_lower_bound) {
            new_lower_bound = cd->diff_children[i].lower_bound;
         }
      }
      if (new_lower_bound > cd->lower_bound) {
         if (! cd->parent) {/* i.e. cd == root_cd */
            printf("Lower bound increased from %d to %d.\n",
                   cd->lower_bound,new_lower_bound);
         }
         cd->lower_bound = new_lower_bound;
         cd = cd->parent;
      } else {
         cd = (colordata*) NULL;
      }
   }   
   return rval;
}


static int remove_finished_subtree(colordata* child)
{
   int rval = 0;
   int i;
   colordata* cd = (colordata*) child;
   int all_same_finished = 1;
   int all_diff_finished = 1;

   while (cd) {
      for (i = 0; i < cd->nsame; ++i) {
         if (cd->same_children[i].status < finished) {
            all_same_finished = 0;break;
         }
      }
      if (cd->nsame && all_same_finished) {
         rval = collect_same_children(cd);
         COLORcheck_rval(rval,"Failed in collect_same_children");
         for (i = 0;  i < cd->nsame; ++i) {
            free_colordata(cd->same_children + i);
         }
         free(cd->same_children);
         cd->same_children = (colordata*) NULL;
         cd->nsame = 0;
      }


      for (i = 0; i < cd->ndiff; ++i) {
         if (cd->diff_children[i].status < finished) {
            all_diff_finished = 0;break;
         }
      }
      if (cd->ndiff && all_diff_finished) {
         rval = collect_diff_children(cd);
         COLORcheck_rval(rval,"Failed in collect_diff_children");
         for (i = 0;  i < cd->ndiff; ++i) {
            free_colordata(cd->diff_children + i);
         }
         free(cd->diff_children);
         cd->diff_children = (colordata*) NULL;
         cd->ndiff = 0;
      }

/*       if (cd->diff_children && cd->diff_children->status == finished) { */
/*          rval = collect_diff_children(cd); */
/*          COLORcheck_rval(rval,"Failed in collect_diff_children"); */
/*          free_colordata(cd->diff_children); */
/*          free(cd->diff_children); */
/*          cd->diff_children = (colordata*) NULL; */
/*       } */
      if (!cd->same_children && !cd->diff_children) {
         cd->status = finished;
         cd = cd->parent;
      } else {
         cd = (colordata*) NULL;
      }
   }
 CLEANUP:
   return rval;
}



static int insert_into_branching_heap(COLORNWTHeap* heap, colordata* cd, double key_mult)
{
   int rval = 0;
   int dummy_href;
   /** We use the lower bound adjusted by the depth of a node as a
       heap key.  The idea behind this is: in case of equal bounds
       deeper nodes should be preferred to potentially improve the
       upper bound faster.
    */
   COLORNWT heap_key = (COLORNWT) (cd->dbl_est_lower_bound * key_mult) - cd->depth;

   if (COLORdbg_lvl() > 1) {
      printf("Inserting into branching heap with lb %d (%f) and ub %d at depth %d (id = %d):\n",
             cd->lower_bound,cd->dbl_est_lower_bound,cd->upper_bound,
             cd->depth,
             cd->id );
      
      print_graph_operations(cd);
   }

   if (cd->lower_bound < cd->upper_bound) {
      COLORcheck_rval(rval,"Failed in compute_lower_bound (cd);");
      rval = COLORNWTheap_insert(heap,&dummy_href,
                                 heap_key,
                                 (void*) cd);
      COLORcheck_rval(rval, "Failed to COLORNWTheap_insert");
   } else {
      cd->status = finished;
      printf("Skipping  with lb %d (%f) and ub %d at depth %d (id = %d, "
             "opt_track = %d, unprocessed nodes = %d) (insertion).\n",
             cd->lower_bound,cd->dbl_est_lower_bound,cd->upper_bound,
             cd->depth,
             cd->id, cd->opt_track,COLORNWTheap_size(heap) );
   }

   rval = trigger_lb_changes(cd);
   COLORcheck_rval(rval,"Failed in trigger_lb_changes");

 CLEANUP:
   return rval;
}

static int sequential_branching(COLORNWTHeap* br_heap,
                                COLORNWT global_upper_bound,
                                double   key_mult)
{
   int rval = 0;
   colordata*    cd;
   printf("ENTERED SEQUENTIAL BRANCHING,\n");

   while ( (cd = (colordata*) COLORNWTheap_min(br_heap) ) ) {
      int i;
      cd->upper_bound = global_upper_bound;
      if (cd->lower_bound < cd->upper_bound) {
         printf("Branching with lb %d (est. %f) and ub %d at depth %d (id = %d, "
                "opt_track = %d, unprocessed nodes = %d).\n",
                cd->lower_bound,cd->dbl_est_lower_bound,cd->upper_bound,
                cd->depth,
                cd->id, cd->opt_track,COLORNWTheap_size(br_heap) );

         /* Create children. If the current LP-solution turns out to be intergal,
            cd->upper_bound might decrease!
          */
         rval = create_branches(cd);
         COLORcheck_rval(rval,"Failed to create_branches");


         for (i = 0; i < cd->nsame; ++i) {
            /* printf("Adding same child to heap\n"); */
            rval = insert_into_branching_heap(br_heap,&(cd->same_children[i]),key_mult);
            COLORcheck_rval(rval,"Failed in insert_into_cb_heap");
         }
         for (i = 0; i < cd->ndiff; ++i) {
/*             printf("Adding diff child to heap\n"); */
            rval = insert_into_branching_heap(br_heap,&(cd->diff_children[i]),key_mult);
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
   return rval;
}

static int locate_colordata(colordata** cd, 
                            branching_joblist** joblist,
                            int cd_id)
{
   int rval = 0;
   branching_joblist*  cur_job = (branching_joblist*) NULL;
   branching_joblist* prev_job = (branching_joblist*) NULL;


   cur_job = *joblist;

   while (cur_job && cur_job->cd->id != cd_id) {
      prev_job = cur_job;
      cur_job  = cur_job->next;
   }
   
   if (!cur_job) {
      fprintf(stderr,"Job %d not found!!!\n",cd_id);
      rval = 1; goto CLEANUP;
   }

   *cd = cur_job->cd;
   if (prev_job) {
      prev_job->next = cur_job->next;
   } else {
      *joblist = cur_job->next;
   }

   COLOR_FREE(cur_job,branching_joblist);

 CLEANUP:
   return rval;
}

static int prepend_to_joblist(branching_joblist** joblist,colordata* cd)
{
   int rval = 0;

   branching_joblist* new_job = COLOR_SAFE_MALLOC(1,branching_joblist);
   COLORcheck_NULL(new_job,"Failed to allocate new_job");

   new_job->cd   = cd;
   new_job->next = *joblist;

   *joblist = new_job;

 CLEANUP:
   return rval;
}


static int parallel_branching(COLORNWTHeap* br_heap,
                              double*       child_cputimes,
                              COLORNWT      global_upper_bound,
                              double        key_mult)
{
   int rval = 0;
   colordata*   cd    = (colordata*)    NULL;
   COLOR_SPORT* lport = (COLOR_SPORT *) NULL;
   COLOR_SFILE* s;
   char request, grunt[1024];
   int  i;
   int npending = 0;
   branching_joblist* joblist = (branching_joblist*) NULL;

   *child_cputimes = 0;

   printf("ENTERED PARALLEL BRANCHING,\n");

   lport = COLORsafe_snet_listen (COLOR_BOSS_PORT);
   if (lport == (COLOR_SPORT *) NULL) {
      fprintf (stderr, "COLORsafe_snet_listen failed\n");
      rval = 1; goto CLEANUP;
   }

   do {
      s = COLORsafe_snet_receive (lport);
      if (!s) {
         fprintf (stderr, "COLORsafe_snet_receive failed, ignoring\n");
         continue;
      }

      if (COLORsafe_sread_string (s, grunt, 1023)) {
         fprintf (stderr, "COLORsafe_sread_char string, abort con\n");
         COLORsafe_sclose (s);
         continue;
      }

      if (COLORsafe_sread_char (s, &request)) {
         fprintf (stderr, "COLORsafe_sread_char failed, abort con\n");
         COLORsafe_sclose (s);
         continue;
      }

      switch (request) {
      case COLOR_BOSS_SEND:
         cd = (colordata*) COLORNWTheap_min(br_heap);
         while (cd && cd->lower_bound >= global_upper_bound) {         
            printf("Skipping with lb %d (%f) and ub %d at depth %d (id = %d, "
                   "opt_track = %d, unprocessed nodes = %d).\n",
                   cd->lower_bound,cd->dbl_est_lower_bound,cd->upper_bound,
                   cd->depth,
                   cd->id, cd->opt_track,COLORNWTheap_size(br_heap) );
            cd->status = finished;
            remove_finished_subtree(cd);
            cd = (colordata*) COLORNWTheap_min(br_heap);
         }
         if (cd) {
            cd->upper_bound = global_upper_bound;
            int include_bestcolors = 0;
            printf("Branching with lb %d (%f) and ub %d at depth %d (id = %d, "
                   "opt_track = %d, unprocessed nodes = %d, branching cputime %f).\n",
                   cd->lower_bound,cd->dbl_est_lower_bound,cd->upper_bound,
                   cd->depth,
                   cd->id, cd->opt_track,COLORNWTheap_size(br_heap),*child_cputimes );
            
            /* Create children. If the current LP-solution turns out to be intergal,
               cd->upper_bound might decrease!
            */
            
            rval = COLORsafe_swrite_char (s, COLOR_BOSS_YES);
            COLORcheck_rval (rval, "COLORsafe_swrite_int failed (YES)");
            
            cd->status = submitted_for_branching;
            npending++;
            rval = send_colordata(s,cd,include_bestcolors);
            COLORcheck_rval(rval,"Failed to create_branches");

            rval = prepend_to_joblist(&joblist,cd);
            COLORcheck_rval(rval,"Failed to add_to_joblist");
            
         } else {
            if (npending) {
               rval = COLORsafe_swrite_char (s, COLOR_BOSS_NO);
            
            } else {
               rval = COLORsafe_swrite_char (s, COLOR_BOSS_EXIT);
            }
            if (rval) {
               fprintf (stderr, "BOSS_NONE write failed - abort con\n");
            }
         }
         COLORsafe_sclose (s);
         break;
      case COLOR_BOSS_RECEIVE:
         {
            int cd_id = -1;
            int include_bestcolors = 1;
            int adopt_id           = 0;
            double child_cputime;
            if (COLORsafe_sread_int (s, &cd_id)) {
               fprintf (stderr, "COLORsafe_sread_int failed, abort cd_id\n");
               COLORsafe_sclose (s);
               continue;
            }

            if (COLORsafe_sread_double (s, &child_cputime)) {
               fprintf (stderr, "COLORsafe_sread_double failed, abort child_cputime\n");
               COLORsafe_sclose (s);
               continue;
            }
            *child_cputimes += child_cputime;

            //  locate cd pointing to cd_id
            rval = locate_colordata(&cd, &joblist,cd_id);
            COLORcheck_rval(rval, "Failed in locate_colordata");

            

            free_temporary_data(cd);
            rval = receive_colordata (s, cd,adopt_id,include_bestcolors);
            printf("Received data with ndiff = %d.\n",cd -> ndiff);
            npending--;
            if (rval) {
               fprintf (stderr, "receive_result failed - abort connection\n");
               COLORsafe_sclose (s);
               break;
            } else {
               printf ("Received result %d from %s\n", cd_id, grunt);
               fflush (stdout);
               COLORsafe_sclose (s);
               for (i = 0; i < cd->nsame; ++i) {
/*                   printf("Adding same child to heap\n"); */
                  rval = insert_into_branching_heap(br_heap,&(cd->same_children[i]),key_mult);
                  COLORcheck_rval(rval,"Failed in insert_into_cb_heap");
               }
               for (i = 0; i < cd->ndiff; ++i) {
/*                   printf("Adding diff child to heap\n"); */
                  rval = insert_into_branching_heap(br_heap,&(cd->diff_children[i]),key_mult);
                  COLORcheck_rval(rval,"Failed in insert_into_cb_heap");
               }
               assert (cd->lower_bound <= cd->upper_bound);
                  
               if (global_upper_bound > cd->upper_bound) {
                  global_upper_bound = cd->upper_bound;
                  
               }
               remove_finished_subtree(cd);
            }
         }
         break;
      case COLOR_BOSS_EXIT:
         printf ("Shutting down the test boss\n"); fflush (stdout);
         COLORsafe_sclose (s);
         goto CLEANUP;
      default:
         fprintf (stderr, "Invalid request %c\n", request);
      }
   }

   while (npending || cd);

 CLEANUP:
   COLORsafe_snet_unlisten (lport);

   return rval;
}


int compute_coloring(colordata* root_cd,
                     int        parallel)
{
   int           rval = 0;
   COLORNWTHeap* br_heap            = (COLORNWTHeap*) NULL;
   double        key_mult           = (double) (COLORNWT_MAX - 1) / root_cd->ncount;
   COLORNWT      global_upper_bound = COLORNWT_MAX;

   double colheur_rtime;
   double branching_rtime;
   double init_lb_rtime;

   init_lb_rtime = - COLORcpu_time();
   COLORNWTheap_init(&br_heap,root_cd->ncount * root_cd->ncount);

   rval = compute_lower_bound(root_cd);
   COLORcheck_rval(rval,"Failed in compute_lower_bound");

   rval = insert_into_branching_heap(br_heap,root_cd,key_mult);
   COLORcheck_rval(rval,"Failed in insert_into_branching_heap");

   rval = write_mwis_instances(root_cd);
   COLORcheck_rval(rval,"Failed ing write_mwis_instances");

   if (cclasses_outfile) {
      int add_timestamp = 0;
      write_snapshot(root_cd,add_timestamp);
      COLORcheck_rval(rval,"Failed ing write_snapshot");
   }
   init_lb_rtime += COLORcpu_time();

   colheur_rtime = -COLORcpu_time();
   heur_colors_with_stable_sets(root_cd);
   colheur_rtime += COLORcpu_time();

   global_upper_bound = (global_upper_bound > root_cd->upper_bound) ?
      root_cd->upper_bound :global_upper_bound;

   printf("Upper bound heuristic on root node took %f seconds.\n",colheur_rtime);

   branching_rtime = -COLORcpu_time();
   
   if (parallel) {
      double child_cputimes;
         rval = parallel_branching(br_heap,&child_cputimes,
                                global_upper_bound,key_mult);
      COLORcheck_rval(rval,"Failed in parallel_branching");
      branching_rtime += child_cputimes;
   } else {
      sequential_branching(br_heap,global_upper_bound,key_mult);
      COLORcheck_rval(rval,"Failed in sequential_branching");
   }
   branching_rtime += COLORcpu_time();

   printf("Compute_coloring took %f seconds (initial lower bound:%f, heur. "
          "upper bound: %f, branching: %f.\n",
          init_lb_rtime + colheur_rtime + branching_rtime,
          init_lb_rtime,colheur_rtime,branching_rtime);
 CLEANUP:
   COLORNWTheap_free(br_heap); br_heap = (COLORNWTHeap*) NULL;

   return rval;
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

