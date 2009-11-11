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

struct _MWISgrb_env {
   GRBenv*   grb_env;
   GRBmodel* grb_model;
   double*   x_opt;
};

/* the actual implementation of gurobi-based MWIS.*/
static int mwis_optimize_model(MWISgrb_env** env,
                               COLORset** newsets, int* nnewsets, int ncount, 
                               double nodelimit,
                               int ecount, const int elist[], double nweights[]);

static int intercept_grb_cb(GRBmodel *grb_model, void *cbdata, int where, void *usrdata)
{
   int rval = 0;
   
   if (where ==GRB_CB_MIPSOL) {
      double objective;
      
      rval = GRBcbget(cbdata,where,GRB_CB_MIPSOL_OBJBST,(void*) &objective);
      COLORcheck_rval (rval, "GRBcbget OBJBST failed");
      
      if (objective > 1 + int_tolerance) {
         if(COLORdbg_lvl()) {
            printf("Terminating gurobi based on current objective value %f\n",
                   objective);
         }
         GRBterminate(grb_model);
      }
   }
   
 CLEANUP:
   return rval;
}


static int mwis_init_model(MWISgrb_env** env,
                           double nodelimit, int ncount,
                           int ecount, const int elist[], double nweights[])
{
   int rval = 0,i;
   char     *vtype = (char *) NULL;     /* variable types*/

   if (!* env) {
      GRBenv* grb_env = (GRBenv*) NULL;
      GRBmodel* grb_model = (GRBmodel*) NULL;

      *env = (MWISgrb_env*) malloc(sizeof(MWISgrb_env));
      COLORcheck_NULL(*env,"Allocating *env failed.");
      (*env)->grb_env   = (GRBenv*) NULL;
      (*env)->grb_model = (GRBmodel*) NULL;
      (*env)->x_opt = (double*) NULL;
      
      rval = GRBloadenv (&((*env)->grb_env), "mwis_gurobi.log");
      if (rval) {
         fprintf (stderr, "GRBloadenv failed\n");
         goto CLEANUP;
      }

      /* Set to 1 to turn on Gurobi output, 0 to turn off output */
      rval = GRBsetintparam ((*env)->grb_env, GRB_INT_PAR_OUTPUTFLAG, COLORdbg_lvl() ? 1 : 0);
      COLORcheck_rval (rval, "GRBsetintparam OUTPUTFLAG failed");
      grb_env = (*env)->grb_env;

      rval = GRBsetintparam (grb_env,GRB_INT_PAR_THREADS , 1);
      if (rval ) {
         fprintf (stderr, "GRBsetintparam GRB_INT_PAR_THREADS failed: %s\n",
                  GRBgeterrormsg(grb_env));
         goto CLEANUP;
      }

      /* Clique cuts should be helpful here, though I didn't see much of
         a difference.
      */
      rval = GRBsetintparam (grb_env,GRB_INT_PAR_CLIQUECUTS , 2);
      if (rval ) {
         fprintf (stderr, "GRBsetintparam GRB_INT_PAR_CLIQUECUTS failed: %s\n",
                  GRBgeterrormsg(grb_env));
         goto CLEANUP;
      }
      
      rval = GRBsetdblparam (grb_env,GRB_DBL_PAR_NODELIMIT , nodelimit);
      if (rval ) {
         fprintf (stderr, "GRBsetdblparam GRB_DBL_PAR_NODELIMIT failed: %s\n",
                  GRBgeterrormsg(grb_env));
         goto CLEANUP;
      }
      
      vtype = (char *) malloc (ncount * sizeof (char));
      if (!vtype) {
         fprintf (stderr, "out of memory for vtype\n");
         rval = 1;  goto CLEANUP;
      }
      for (i = 0; i < ncount; i++) vtype[i] = GRB_BINARY;
      
      (*env)->x_opt = (double *) malloc (ncount * sizeof (double));
      if (!(*env)->x_opt) {
         fprintf (stderr, "out of memory for vtype\n");
         rval = 1;  goto CLEANUP;
      }

      rval = GRBnewmodel (grb_env, &((*env)->grb_model), 
                          "mwisme", ncount, (double*) nweights,
                          (double *) NULL, (double *) NULL, vtype, (char**) NULL);
      if (rval) {
         fprintf (stderr, "GRBnewmodel failed: %s\n",
                  GRBgeterrormsg(grb_env));
         goto CLEANUP;
      }
      grb_model = (*env)->grb_model;
      
      rval  = GRBsetcallbackfunc(grb_model, intercept_grb_cb, NULL);
      COLORcheck_rval (rval, "GRBsetcallbackfunc failed");


      /* We are dealing with a maximization  problem. */
      rval = GRBsetintattr(grb_model,GRB_INT_ATTR_MODELSENSE, -1);
      if (rval) {
         fprintf (stderr, "GRBsetintattr: %s\n",
                  GRBgeterrormsg(grb_env));
         goto CLEANUP;
      }
   
      for (i = 0; i < ecount; i++) {
         int numnnz     = 2;
         int v[2];
         double coef[2] = {1.0, 1.0};
         
         v[0]  = elist[2*i];
         v[1]  = elist[2*i +1];
         
         rval = GRBaddconstr (grb_model, numnnz,v,coef,
                              GRB_LESS_EQUAL, 1.0, NULL);
         if (rval) {
            fprintf (stderr, "MWIS GRBaddconstr failed: %s (i: %d u: %d v: %d)\n",
                     GRBgeterrormsg(grb_env),i,v[0],v[1]);
            goto CLEANUP;
         }
      }
      
      rval = GRBupdatemodel (grb_model);
      if (rval) {
         fprintf (stderr, "GRPupdatemodel failed: %s\n",
                  GRBgeterrormsg(grb_env));
         goto CLEANUP;
      }

#ifdef WRITE_MODEL   
      rval = GRBwrite (grb_model, "mwispre.lp");
      if (rval) {
         fprintf (stderr, "GRPwrite failed: %s\n",
                  GRBgeterrormsg(grb_env));
         goto CLEANUP;
      }
#endif 
   }
 CLEANUP:
   if (vtype) free(vtype);
   
   if(rval) {
      int frval = COLORstable_free_grb_env(env);
      COLORcheck_rval (frval, "COLORstable_free_grb_env failed");
   }
   return rval;
}


/** Solve the maximum weighted independent set problem given
    by ncount, ecount, elist, nweights, with gurobi and store newly
    generated solution(s) in newsets, nnewsets.
   
 */
int COLORstable_gurobi(MWISgrb_env** env,
                       COLORset** newsets, int* nnewsets, int ncount, 
    int ecount, const int elist[], double nweights[])
{
   int    rval;
   double nodelimit = DBL_MAX;

   assert(*newsets == (COLORset*) NULL);
   
   rval = mwis_optimize_model(env,
                              newsets,nnewsets,
                              ncount,nodelimit,ecount,elist,
                              nweights);
   if (rval) {
      fprintf (stderr, "mwis_optimize_model  failed.\n");
      goto CLEANUP;
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


static int mwis_optimize_model(MWISgrb_env** env,
                               COLORset** newsets, int* nnewsets, int ncount, 
                               double nodelimit,
                               int ecount, const int elist[], double nweights[])
{
   int       rval = 0,i,j;
   double    objective;
   
   if (!*env) {
      rval = mwis_init_model(env,nodelimit,ncount,ecount,elist,nweights);
      COLORcheck_rval (rval, "mwis_init_model failed");
   } else { /* Set to new objective. */
      rval = GRBsetdblattrarray((*env)->grb_model,GRB_DBL_ATTR_OBJ,
                                0,ncount,nweights);
   }
   

   

   rval = GRBoptimize((*env)->grb_model);
   if (rval) {
      fprintf (stderr, "GRBoptimize failed: %s\n",
               GRBgeterrormsg((*env)->grb_env));
      goto CLEANUP;
   }
   
   rval = GRBgetdblattr((*env)->grb_model, GRB_DBL_ATTR_OBJVAL,
                        &objective);
   if (rval) {
      fprintf (stderr, "GRBgetdblattr GRB_DBL_ATTR_OBJVAL failed: %s\n",
               GRBgeterrormsg((*env)->grb_env));
      goto CLEANUP;
   }
   
   if (objective > 1.0 + int_tolerance) {
      COLORset* newset = (COLORset *) NULL;

      /* Retrieve variable values.*/
      rval = GRBgetdblattrarray((*env)->grb_model, GRB_DBL_ATTR_X,
                                0, ncount,(*env)->x_opt);
      if (rval) {
         fprintf (stderr, "GRBgetdblattrarray failed: %s\n",
                  GRBgeterrormsg((*env)->grb_env));
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
         if ( value_is_one((*env)->x_opt[i]) ) {
            newset->count++;
         } else if (fabs (!value_is_zero((*env)->x_opt[i]))) {
            /* we better should check wether we found a desired solution,
               i.e. x_opt defines an independent set with "nweights*x_opt > 1.0".
            */
            fprintf (stderr, "Found non-binary solution: var: %d x: %g\n",
                     i, (*env)->x_opt[i]);
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
         if ( value_is_one((*env)->x_opt[i]) ) {
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

   if (rval) {
      int frval = COLORstable_free_grb_env(env);
      if(frval) printf("COLORstable_free_grb_env failed!\n");
      COLORfree_sets (newsets,nnewsets);
   }                                            
   return rval;
}


int COLORstable_free_grb_env(MWISgrb_env** env)
{
   int rval = 0;
   if (env) {
      if (*env) {
         if ((*env)->x_opt) free ((*env)->x_opt);
         if ((*env)->grb_model) {
            rval = GRBfreemodel((*env)->grb_model);
            if(rval) printf("GRBfreemodel failed.\n");
         }
         if((*env)->grb_env) GRBfreeenv((*env)->grb_env);
      
         free(*env);
         *env = (MWISgrb_env*) NULL;
      }
   }
   return rval;
}


int COLORstable_write_mps(const char*  filename,
                          int ncount, int ecount, const int elist[], double nweights[])
{
   int rval = 0;
   MWISgrb_env* env = (MWISgrb_env*) NULL;
   rval = mwis_init_model(&env,1,ncount,ecount,elist,nweights);
   COLORcheck_rval(rval,"Failed in mwis_init_model");

   
   rval = GRBwrite (env->grb_model, filename);
   if (rval) {
      fprintf (stderr, "GRPwrite failed: %s\n",
               GRBgeterrormsg(env->grb_env));
      goto CLEANUP;
   }

 CLEANUP:
   return rval;
}
