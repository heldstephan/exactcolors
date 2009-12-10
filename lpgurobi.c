#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lp.h"
#include "color.h"


const double int_tolerance = 0.00001;

#define COLORcheck_rval_grb(rval,msg,env) {                            \
      if ((rval)) {                                                    \
         fprintf (stderr, "%s at %s, line %d: %s\n",                   \
                  (msg), __FILE__, __LINE__,GRBgeterrormsg (env));     \
         goto CLEANUP;                                                 \
    }                                                                  \
}


int COLORlp_init (COLORlp **p, const char *name)
{
    int rval = 0;

    (*p) = (COLORlp *) COLOR_SAFE_MALLOC (1,COLORlp);
    COLORcheck_NULL (*p, "out of memory for lp");

    (*p)->env = (GRBenv *) NULL;
    (*p)->model = (GRBmodel *) NULL;

    rval = GRBloadenv (&((*p)->env), NULL);
    COLORcheck_rval (rval, "GRBloadenv failed");

    /* Set to 1 to turn on Gurobi output, 0 to turn off output */
    rval = GRBsetintparam ((*p)->env, GRB_INT_PAR_OUTPUTFLAG, COLORdbg_lvl() ? 1 : 0);
    COLORcheck_rval_grb (rval, "GRBsetintparam OUTPUTFLAG failed",(*p)->env);

    rval = GRBsetintparam ((*p)->env, GRB_INT_PAR_THREADS , 1);
    COLORcheck_rval_grb (rval, "GRBsetintparam THREADS failed",(*p)->env);

    rval = GRBnewmodel ((*p)->env, &((*p)->model), name, 0, (double *) NULL,
                    (double *) NULL, (double *) NULL, (char *) NULL, NULL);
    COLORcheck_rval_grb (rval, "GRBnewmodel failed",(*p)->env);

CLEANUP:
    return rval;
}

void COLORlp_free (COLORlp **p)
{
    if (*p) {
        if ((*p)->model) GRBfreemodel ((*p)->model);
        if ((*p)->env) GRBfreeenv ((*p)->env);
        free (*p);
        *p = (COLORlp *) NULL;
    }
}

int COLORlp_optimize (COLORlp *p)
{
    int rval = 0;

    rval = GRBoptimize (p->model);
    COLORcheck_rval_grb (rval, "GRBoptimize failed",p->env);

CLEANUP:
    return rval;
}

int COLORlp_objval (COLORlp *p, double *obj)
{
    int rval = 0;

    rval = GRBgetdblattr (p->model, GRB_DBL_ATTR_OBJVAL, obj);
    COLORcheck_rval_grb (rval, "GRBgetdblattr OBJVAL failed",p->env);

CLEANUP:
    return rval;

}

int COLORlp_change_objective(COLORlp *p, int start, int len, double* values)
{
   int rval = 0;
   rval = GRBsetdblattrarray(p->model,GRB_DBL_ATTR_OBJ,
                             start,len,values);
   COLORcheck_rval_grb(rval,"Failed in GRBsetdblattrarray",p->env);

   rval = GRBupdatemodel (p->model);
   COLORcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

 CLEANUP:
   return rval;
}

int COLORlp_addrow (COLORlp *p, int nzcount, int *cind, double *cval, 
       char sense, double rhs, char *name) 
{
    int rval = 0;

    rval = GRBaddconstr (p->model, nzcount, cind, cval, sense, rhs, name); 
    COLORcheck_rval_grb (rval, "GRBaddconstr failed", p->env);
    rval = GRBupdatemodel (p->model);
    COLORcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

CLEANUP:

    return rval;
}

int COLORlp_addcol (COLORlp *p, int nzcount, int *cind, double *cval, 
       double obj, double lb, double ub, char vartype, char *name) 
{
    int rval = 0;

    rval = GRBaddvar (p->model, nzcount, cind, cval, obj, lb, ub, vartype,
                      name); 
    COLORcheck_rval_grb (rval, "GRBaddvar failed", p->env);
    rval = GRBupdatemodel (p->model);
    COLORcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

CLEANUP:

    return rval;
}

int COLORlp_deletecol (COLORlp *p, int cind)
{
   int rval = 0;

   rval = GRBdelvars(p->model,1,&cind);
   COLORcheck_rval_grb (rval, "GRBdelvars failed", p->env);
   rval = GRBupdatemodel (p->model);
   COLORcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);
 CLEANUP:
   
   return rval;
}


int COLORlp_pi (COLORlp *p, double *pi)
{
    int rval = 0;
    int nrows;
    int solstat;

    rval = GRBgetintattr(p->model, GRB_INT_ATTR_STATUS, &solstat);
    COLORcheck_rval_grb (rval, "GRBgetintattr GRB_INT_ATTR_STATUS failed", 
                         p->env);
    if (solstat == GRB_INFEASIBLE) {
        fprintf (stderr, "Problem is infeasible\n");
        rval = 1;  goto CLEANUP;
    }

    rval = GRBgetintattr(p->model, GRB_INT_ATTR_NUMCONSTRS, &nrows);
    if (nrows == 0) {
        fprintf (stderr, "No rows in LP\n");
        rval = 1;  goto CLEANUP;
    }

    rval = GRBgetdblattrarray(p->model, GRB_DBL_ATTR_PI, 0, nrows, pi);
    COLORcheck_rval_grb (rval, "GRBgetdblattrarray GRB_DBL_ATTR_PI failed", 
                         p->env);

CLEANUP:
    return rval;
}

int COLORlp_x (COLORlp *p, double *x)
{
    int rval = 0;
    int ncols;
    int solstat;

    rval = GRBgetintattr(p->model, GRB_INT_ATTR_STATUS, &solstat);
    COLORcheck_rval_grb (rval, "GRBgetintattr GRB_INT_ATTR_STATUS failed", 
                         p->env);
    if (solstat == GRB_INFEASIBLE) {
        fprintf (stderr, "Problem is infeasible\n");
        rval = 1;  goto CLEANUP;
    }

    rval = GRBgetintattr(p->model, GRB_INT_ATTR_NUMVARS, &ncols);
    if (ncols == 0) {
        fprintf (stderr, "No columns in LP\n");
        rval = 1;  goto CLEANUP;
    }

    rval = GRBgetdblattrarray(p->model, GRB_DBL_ATTR_X, 0, ncols, x);
    COLORcheck_rval_grb (rval, "GRBgetdblattrarray GRB_DBL_ATTR_X failed", 
                         p->env);

CLEANUP:
    return rval;
}

int COLORlp_set_all_coltypes (COLORlp *p, char sense)
{
   int nvars,i;
   int rval= GRBgetintattr(p->model,GRB_INT_ATTR_NUMVARS,&nvars);
   COLORcheck_rval_grb (rval, "GRBgetintattr GRB_INT_ATTR_NUMVARS failed",
                        p->env);
   
   for (i = 0; i < nvars; i++) {
      rval = GRBsetcharattrelement(p->model,GRB_CHAR_ATTR_VTYPE,i,sense);
      COLORcheck_rval_grb (rval, "GRBsetintattrelement GRB_CHAR_ATTR_VTYPE failed",
                           p->env);
   }

   rval = GRBupdatemodel (p->model);
   COLORcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

 CLEANUP:
   return rval;
}

int COLORlp_objective_sense (COLORlp *p, int sense)
{
    int rval = 0;

    /* Min = 1   Max = -1 */
    rval = GRBsetintattr (p->model, GRB_INT_ATTR_MODELSENSE, sense);
    COLORcheck_rval_grb (rval, "GRBsetintattr failed", p->env);

CLEANUP:

    return rval;
}

int COLORlp_setbound (COLORlp *p, int col, char lower_or_upper, double bnd)
{
    int rval = 0;

    if (lower_or_upper == 'L') {
        rval = GRBsetdblattrelement (p->model, GRB_DBL_ATTR_LB, col, bnd);
    } else {
        rval = GRBsetdblattrelement (p->model, GRB_DBL_ATTR_UB, col, bnd);
    }
    COLORcheck_rval (rval, "GRBsetdblattr LB or UB failed");

    rval = GRBupdatemodel (p->model);
    COLORcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

CLEANUP:

    return rval;
}

int COLORlp_setnodelimit (COLORlp *p, int mip_node_limit)
{
   int rval = GRBsetdblparam (GRBgetenv(p->model), GRB_DBL_PAR_NODELIMIT, mip_node_limit);
   COLORcheck_rval_grb (rval, "GRBsetdblparam NODELIMIT failed",p->env);
 CLEANUP:
   return rval;
}

int COLORlp_write (COLORlp *p, const char *fname)
{
    int rval = 0;

    rval = GRBwrite (p->model, fname);
    COLORcheck_rval_grb (rval, "GRBwrite failed", p->env);

CLEANUP:

    return rval;
}

void COLORlp_free_warmstart (COLORlp_warmstart **w)
{
    if (*w != (COLORlp_warmstart *) NULL) {
        COLOR_IFFREE ((*w)->cstat, int);
        COLOR_IFFREE ((*w)->rstat, int);
        COLOR_IFFREE ((*w)->dnorm, double);
        COLOR_FREE (*w, COLORlp_warmstart);
    }
}

void COLORlp_printerrorcode (int c)
{
    switch (c) {
    case GRB_ERROR_OUT_OF_MEMORY:
        printf ("Available memory was exhausted\n");
        break;
    case GRB_ERROR_NULL_ARGUMENT:
        printf ("NULL input value provided for a required argument\n");
        break;
    case GRB_ERROR_INVALID_ARGUMENT:
        printf ("An invalid value was provided for a routine argument\n");
        break;
    case GRB_ERROR_UNKNOWN_ATTRIBUTE:
        printf ("Tried to query or set an unknown attribute\n");
        break;
    case GRB_ERROR_DATA_NOT_AVAILABLE:
        printf ("Attempted to query or set an attribute that could\n");
        printf ("not be accessed at that time\n");
        break;
    case GRB_ERROR_INDEX_OUT_OF_RANGE:
        printf ("Tried to query or set an attribute, but one or more\n");
        printf ("of the provided indices (e.g., constraint index, variable \n");
        printf ("index) was outside the range of valid values\n");
        break;
    case GRB_ERROR_UNKNOWN_PARAMETER:
        printf ("Tried to query or set an unknown parameter\n");
        break;
    case GRB_ERROR_VALUE_OUT_OF_RANGE:
        printf ("Tried to set a parameter to a value that is outside\n");
        printf ("the parameter's valid range\n");
        break;
    case GRB_ERROR_NO_LICENSE:
        printf ("Failed to obtain a valid license\n");
        break;
    case GRB_ERROR_SIZE_LIMIT_EXCEEDED:
        printf ("Attempted to solve a model that is larger than the\n");
        printf ("limit for a demo license\n");
        break;
    case GRB_ERROR_CALLBACK:
        printf ("Problem in callback\n");
        break;
    case GRB_ERROR_FILE_READ:
        printf ("Failed to read the requested file\n");
        break;
    case GRB_ERROR_FILE_WRITE:
        printf ("Failed to write the requested file\n");
        break;
    case GRB_ERROR_NUMERIC:
        printf ("Numerical error during requested operation\n");
        break;
    case GRB_ERROR_IIS_NOT_INFEASIBLE:
        printf ("Attempted to perform infeasibility analysis on a\n");
        printf ("feasible model\n"); 
        break;
    case GRB_ERROR_NOT_FOR_MIP:
        printf ("Requested operation not valid for a MIP model\n");
        break;
    case GRB_ERROR_OPTIMIZATION_IN_PROGRESS:
        printf ("Tried to query or modify a model while optimization\n");
        printf ("was in progress\n");
        break;
     default:
        printf ("Unknown error code: %d\n", c);
    }
    fflush (stdout);
}

double COLORlp_int_tolerance ()
{
   return int_tolerance;
}
