#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <gurobi_c.h>
#include "color.h"
#include "mwis.h"


/** Maximum-weight stable-set problem via MIP code (Gurobi) **/
/** Author:  Stephan Held, 091102                           **/

const double int_tolerance = 0.0000001;


/* the actual implementation of gurobi-based MWIS.*/
static int mwis_optimize_model(COLORset** newsets, int* nnewsets, int ncount, 
                               double nodelimit,
                               int ecount, const int elist[], double nweights[]);


/** Solve the maximum weighted independent set problem given
    by ncount, ecount, elist, nweights, with gurobi and store newly
    generated solution(s) in newsets, nnewsets.
   
 */
int COLORstable_gurobi(COLORset** newsets, int* nnewsets, int ncount, 
    int ecount, const int elist[], double nweights[])
{
   int    rval;
   double nodelimit = 1.0;/* B&B nodelimit.*/

   assert(*newsets == (COLORset*) NULL);
   
   /* First try w/o branching.*/
   rval = mwis_optimize_model(newsets,nnewsets,ncount,nodelimit,ecount,elist,
                             nweights);
   if (rval ) {
      fprintf (stderr, "mwis_optimize_model 1 failed.\n");
      goto CLEANUP;
   }
   
   if (! *newsets) {
      /* no success without branching => try with  branching:*/
      nodelimit = DBL_MAX;
      printf("Trying 2nd MWIS.\n");
      rval = mwis_optimize_model(newsets,nnewsets,ncount,nodelimit,ecount,elist,
                                 nweights);
      if (rval) {
         fprintf (stderr, "mwis_optimize_model 2 failed.\n");
         goto CLEANUP;
      }
   }
 CLEANUP:
   return rval;
}


static int value_is_one(double v)
{
   if (fabs(v- 1.0) < int_tolerance ) {
      return 1;
   }
   return 0;
}

static int value_is_zero(double v)
{
   if (fabs (v) < int_tolerance) {
      return 1;
   }
   return 0;
}


static int mwis_optimize_model(COLORset** newsets, int* nnewsets, int ncount, 
                               double nodelimit,
                               int ecount, const int elist[], double nweights[])
{
   int rval,gval,i,j;

   GRBmodel *model = (GRBmodel *) NULL;
   char     *vtype = (char *) NULL;     /* variable types*/
   double   *x_opt = (double *) NULL;   /* solution vector*/
   double    objective;

   GRBenv* env;

   rval = GRBloadenv (&env, "mwis_gurobi.log");
   if (rval) {
      fprintf (stderr, "GRBloadenv failed\n");
      goto CLEANUP;
   }

   /* Set to 1 to turn on Gurobi output, 0 to turn off output */
   rval = GRBsetintparam (env, GRB_INT_PAR_OUTPUTFLAG, COLORdbg_lvl() ? 1 : 0);
   COLORcheck_rval (rval, "GRBsetintparam OUTPUTFLAG failed");

#ifdef SINGLE_THREAD   
   rval = GRBsetintparam (env,GRB_INT_PAR_THREADS , 1);
   if (rval ) {
      fprintf (stderr, "GRBsetintparam GRB_INT_PAR_THREADS failed: %s\n",
               GRBgeterrormsg(env));
      goto CLEANUP;
   }
#endif

   /* Clique cuts should be helpful here, though I didn't see much of
      a difference.
   */
   rval = GRBsetintparam (env,GRB_INT_PAR_CLIQUECUTS , 2);
   if (rval ) {
      fprintf (stderr, "GRBsetintparam GRB_INT_PAR_CLIQUECUTS failed: %s\n",
               GRBgeterrormsg(env));
      goto CLEANUP;
   }

   rval = GRBsetdblparam (env,GRB_DBL_PAR_NODELIMIT , nodelimit);
   if (rval ) {
      fprintf (stderr, "GRBsetdblparam GRB_DBL_PAR_NODELIMIT failed: %s\n",
               GRBgeterrormsg(env));
      goto CLEANUP;
   }

   vtype = (char *) malloc (ncount * sizeof (char));
   if (!vtype) {
      fprintf (stderr, "out of memory for vtype\n");
      rval = 1;  goto CLEANUP;
   }
   for (i = 0; i < ncount; i++) vtype[i] = GRB_BINARY;

   x_opt = (double *) malloc (ncount * sizeof (double));
   if (!x_opt) {
      fprintf (stderr, "out of memory for vtype\n");
      rval = 1;  goto CLEANUP;
   }
   
   rval = GRBnewmodel (env, &model, "mwisme", ncount, (double*) nweights,
                       (double *) NULL, (double *) NULL, vtype, (char**) NULL);
   if (rval) {
      fprintf (stderr, "GRBnewmodel failed: %s\n",
               GRBgeterrormsg(env));
      goto CLEANUP;
   }

   /* We are dealing with a maximization  problem. */
   rval = GRBsetintattr(model,GRB_INT_ATTR_MODELSENSE, -1);
   if (rval) {
      fprintf (stderr, "GRBsetintattr: %s\n",
               GRBgeterrormsg(env));
      goto CLEANUP;
   }
   
   for (i = 0; i < ecount; i++) {
      int numnnz     = 2;
      int v[2];
      double coef[2] = {1.0, 1.0};
      
      v[0]  = elist[2*i];
      v[1]  = elist[2*i +1];
      
      rval = GRBaddconstr (model, numnnz,v,coef,
                           GRB_LESS_EQUAL, 1.0, NULL);
      if (rval) {
         fprintf (stderr, "MWIS GRBaddconstr failed: %s (i: %d u: %d v: %d)\n",
                  GRBgeterrormsg(env),i,v[0],v[1]);
         goto CLEANUP;
      }
   }

   rval = GRBupdatemodel (model);
   if (rval) {
      fprintf (stderr, "GRPupdatemodel failed: %s\n",
               GRBgeterrormsg(env));
      goto CLEANUP;
   }

#ifdef WRITE_MODEL   
   rval = GRBwrite (model, "mwispre.lp");
   if (rval) {
      fprintf (stderr, "GRPwrite failed: %s\n",
               GRBgeterrormsg(env));
      goto CLEANUP;
   }
#endif 

   rval = GRBoptimize(model);
   if (rval) {
      fprintf (stderr, "GRBoptimize failed: %s\n",
               GRBgeterrormsg(env));
      goto CLEANUP;
   }
   
   rval = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL,
                        &objective);
   if (rval) {
      fprintf (stderr, "GRBgetdblattr GRB_DBL_ATTR_OBJVAL failed: %s\n",
               GRBgeterrormsg(env));
      goto CLEANUP;
   }
   
   if (objective > 1.0 + int_tolerance) {
      COLORset* newset = (COLORset *) NULL;

      /* Retrieve variable values.*/
      rval = GRBgetdblattrarray(model, GRB_DBL_ATTR_X,
                                0, ncount,x_opt);
      if (rval) {
         fprintf (stderr, "GRBgetdblattrarray failed: %s\n",
                  GRBgeterrormsg(env));
         goto CLEANUP;
      }

      /* Currently we only retrieve a single set.*/
      *nnewsets = 1;    
      *newsets = (COLORset *) malloc(sizeof(COLORset));     

      if (! *newsets) {
         fprintf (stderr, "out of memory for newsets\n");
         rval = 1;  goto CLEANUP;
      }

      /* Firstly, count # of 1-values.*/
      newset = &((*newsets)[0]);
      newset->count = 0;
      for (i = 0; i < ncount;++i) {
         if ( value_is_one(x_opt[i]) ) {
            newset->count++;
         } else if (fabs (!value_is_zero(x_opt[i]))) {
            /* we better should check wether we found a desired solution,
               i.e. x_opt defines an independent set with "nweights*x_opt > 1.0".
            */
            fprintf (stderr, "Found non-binary solution: var: %d x: %g\n",
                     i, x_opt[i]);
            rval = 1; goto CLEANUP;
         }
      }
   
      /* Secondly, generate new independent set.*/
      newset->members = (int *) malloc(newset->count * sizeof(int));
      if (!newset->members) {
         fprintf (stderr, "out of memory for newset.members\n");
         rval = 1;  goto CLEANUP;
      }
      printf("NEW SET ");
      for (i = 0, j = 0; i < ncount;++i) {
         if ( value_is_one(x_opt[i]) ) {
            newset->members[j++] = i;
            printf(" %d",i);
         }
      }
      printf("\n");
   } else if (objective >  1.0 - int_tolerance) {
      fprintf(stderr,"WARNING: MWIS is hardly decidable with objective %g.\n",
              objective);
      goto CLEANUP;
   }
   

 CLEANUP:
   if (model) {
      gval = GRBfreemodel(model);
      if (gval) {
         fprintf (stderr, "MWIS GRBfreemodel failed: %s\n",
                  GRBgeterrormsg(env));
      }
   }
   if (vtype) {free (vtype);}
   if (x_opt) {free (x_opt);}
   if (rval) {
      COLORfree_sets (newsets,nnewsets);
   }                                            

   GRBfreeenv(env);

   return rval;

}
