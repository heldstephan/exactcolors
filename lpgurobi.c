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
#include "lp.h"
#include "color.h"
#include <gurobi_c.h>

struct COLORlp {
    GRBenv *env;
    GRBmodel *model;
};

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
    rval = GRBsetintparam ((*p)->env, GRB_INT_PAR_OUTPUTFLAG, (COLORdbg_lvl() > 0) ? 1 : 0);
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
    int status;

    rval = GRBoptimize (p->model);
    COLORcheck_rval_grb (rval, "GRBoptimize failed",p->env);
    
    rval = GRBgetintattr(p->model,GRB_INT_ATTR_STATUS,&status);
    COLORcheck_rval_grb (rval, "GRBgetintattr failed",p->env);

    if (status != GRB_OPTIMAL) {
       printf("Failed to solve model to optimality. status = ");
       switch (status) {
       case GRB_LOADED:
          printf("GRB_LOADED ");
          rval = 1;break;
       case GRB_INFEASIBLE:
          printf("GRB_INFEASIBLE ");
          rval = GRBcomputeIIS(p->model);
          COLORcheck_rval_grb (rval, "GRBcomputeIIS failed",p->env);
          rval = GRBwrite(p->model,"grbinfeas_debug.lp");
          COLORcheck_rval_grb (rval, "GRBwrite lp failed",p->env);
          
          rval = 1;break;
       case GRB_INF_OR_UNBD:
          printf("GRB_INF_OR_UNBD ");
          rval = 1;break;
       case GRB_UNBOUNDED:
          printf("GRB_UNBOUNDED ");
          rval = 1;break;
       default:
          printf("%d",status);
       }
       printf("\n");

       if(rval) { goto CLEANUP; }
    }

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
    char isense;

    switch (sense) {
    case COLORlp_EQUAL: 
        isense = GRB_EQUAL; break;
    case COLORlp_LESS_EQUAL:    
        isense = GRB_LESS_EQUAL; break;
    case COLORlp_GREATER_EQUAL:    
        isense = GRB_GREATER_EQUAL; break;
    default: 
        fprintf (stderr, "unknown variable sense: %c\n", sense);
        rval = 1;  goto CLEANUP;
    }

    rval = GRBaddconstr (p->model, nzcount, cind, cval, isense, rhs, name); 
    COLORcheck_rval_grb (rval, "GRBaddconstr failed", p->env);
    rval = GRBupdatemodel (p->model);
    COLORcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

CLEANUP:

    return rval;
}

int COLORlp_addcol (COLORlp *p, int nzcount, int *cind, double *cval, 
       double obj, double lb, double ub, char sense, char *name) 
{
    int rval = 0;
    char isense;

    switch (sense) {
    case COLORlp_CONTINUOUS: 
        isense = GRB_CONTINUOUS; break;
    case COLORlp_BINARY:    
        isense = GRB_BINARY; break;
    case COLORlp_INTEGER:    
        isense = GRB_INTEGER; break;
    default: 
        fprintf (stderr, "unknown variable sense: %c\n", sense);
        rval = 1;  goto CLEANUP;
    }

    rval = GRBaddvar (p->model, nzcount, cind, cval, obj, lb, ub, isense,
                      name); 
    COLORcheck_rval_grb (rval, "GRBaddvar failed", p->env);
    rval = GRBupdatemodel (p->model);
    COLORcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

CLEANUP:

    return rval;
}

int COLORlp_deletecols (COLORlp *p, int first_cind, int last_cind)
{
   int rval = 0;
   int* dellist = (int*) NULL;
   int numdel  = last_cind - first_cind + 1;
   int i;
   
   dellist = COLOR_SAFE_MALLOC(numdel,int);
   COLORcheck_NULL(dellist, "Failed to allocate dellist");
   
   for (i = 0; i < numdel; ++i) {
      dellist[i] = first_cind + i;
   }

   rval = GRBdelvars(p->model,numdel,dellist);
   COLORcheck_rval_grb (rval, "GRBdelvars failed", p->env);
   rval = GRBupdatemodel (p->model);
   COLORcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);
 CLEANUP:
   if(dellist) {free (dellist);}

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
    COLORcheck_rval_grb (rval, "GRBgetintattr GRB_INT_ATTR_NUMCONSTRS failed", 
                         p->env);

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
    COLORcheck_rval(rval,"Failed in GRBgetintattr");

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

int COLORlp_basis_cols (COLORlp *p, int *cstat)
{
   int rval = 0;
   int ncols,i;

   rval = GRBgetintattr(p->model, GRB_INT_ATTR_NUMVARS, &ncols);
   COLORcheck_rval(rval,"Failed in GRBgetintattr");

   rval = GRBgetintattrarray(p->model, GRB_INT_ATTR_VBASIS,0,ncols,cstat);
   COLORcheck_rval(rval,"Failed in GRBgetintattrarray");
   
   for (i = 0; i < ncols; ++i) {
      switch (cstat[i]) {
      case 0:
         cstat[i] = COLORlp_BASIC;
         break;
      case -1:
         cstat[i] = COLORlp_LOWER;
         break;
      case -2:
         cstat[i] = COLORlp_UPPER;
         break;
      case -3:
         cstat[i] = COLORlp_FREE;
         break;
      default:
         rval = 1;
         COLORcheck_rval(rval,"ERROR: Received unknown cstat");
      }
   }
 CLEANUP:
   return rval;
}

int COLORlp_set_all_coltypes (COLORlp *p, char sense)
{
   int nvars,i, rval = 0;
   char isense;

   switch (sense) {
   case COLORlp_CONTINUOUS: 
       isense = GRB_CONTINUOUS; break;
   case COLORlp_BINARY:    
       isense = GRB_BINARY; break;
   case COLORlp_INTEGER:    
       isense = GRB_INTEGER; break;
   default: 
       fprintf (stderr, "unknown variable sense: %c\n", sense);
       rval = 1;  goto CLEANUP;
   }

   rval= GRBgetintattr(p->model,GRB_INT_ATTR_NUMVARS,&nvars);
   COLORcheck_rval_grb (rval, "GRBgetintattr GRB_INT_ATTR_NUMVARS failed",
                        p->env);
   
   for (i = 0; i < nvars; i++) {
      rval = GRBsetcharattrelement(p->model,GRB_CHAR_ATTR_VTYPE,i,isense);
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
