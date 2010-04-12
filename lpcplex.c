/****************************************************************************
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
****************************************************************************/

/***  Interface for ILOG CPLEX 12.1  ***/

#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include "lp.h"
#include "color.h"
#include <cplex.h>

struct COLORlp {
    CPXENVptr cplex_env;
    CPXLPptr  cplex_lp;
    int       noptcalls;
};

const double int_tolerance = 0.00001;

int COLORlp_init (COLORlp **p, const char *name)
{
    int rval = 0;


    (*p) = COLOR_SAFE_MALLOC (1, COLORlp);
    if ((*p) == (COLORlp *) NULL) {
        fprintf (stderr, "Out of memory in COLORlp_init\n");
        rval = 1; goto CLEANUP;
    }
    (*p)->noptcalls = 0;
    (*p)->cplex_env = (CPXENVptr) NULL;
    (*p)->cplex_lp = (CPXLPptr) NULL;

    (*p)->cplex_env = CPXopenCPLEX (&rval);
    if (rval) {
        fprintf (stderr, "CPXopenCPLEX failed, return code %d\n", rval);
        goto CLEANUP;
    }

    rval = CPXsetintparam ((*p)->cplex_env, CPX_PARAM_SCRIND, 
                           (COLORdbg_lvl() > 1) ? CPX_ON : CPX_OFF);
    COLORcheck_rval (rval, "CPXsetintparam CPX_PARAM_SCRIND failed");


    rval = CPXsetintparam ((*p)->cplex_env, CPX_PARAM_ADVIND, 1);
    COLORcheck_rval (rval, "CPXsetintparam CPX_PARAM_ADVIND failed");

    rval = CPXsetintparam ((*p)->cplex_env, CPX_PARAM_DPRIIND,
                           CPX_DPRIIND_STEEP);
    COLORcheck_rval (rval, "CPXsetintparam CPX_PARAM_DPRIIND failed");

    rval = CPXsetintparam ((*p)->cplex_env, CPX_PARAM_PPRIIND,
                           CPX_PPRIIND_STEEP);
    COLORcheck_rval (rval, "CPXsetintparam CPX_PARAM_PPRIIND failed");

    /* The following three parameters were set by Bix in TSP code */

    rval = CPXsetdblparam ((*p)->cplex_env, CPX_PARAM_EPPER, 1.0E-6);
    COLORcheck_rval (rval, "CPXsetintparam CPX_PARAM_EPPER failed");

    rval = CPXsetdblparam ((*p)->cplex_env, CPX_PARAM_EPOPT, 1.0E-9);
    COLORcheck_rval (rval, "CPXsetintparam CPX_PARAM_EPOPT failed");

    rval = CPXsetdblparam ((*p)->cplex_env, CPX_PARAM_EPRHS, 1.0E-9);
    COLORcheck_rval (rval, "CPXsetintparam CPX_PARAM_EPRHS failed");


    (*p)->cplex_lp = CPXcreateprob ((*p)->cplex_env, &rval, name);
    if (!(*p)->cplex_lp || rval) {
       fprintf (stderr, "CPXcreateprob failed, return code %d\n", rval);
       goto CLEANUP;
    }

CLEANUP:
    return rval;
}

void COLORlp_free (COLORlp **p)
{
     if (*p) {
        if ((*p)->cplex_env) {
            if ((*p)->cplex_lp) {
                CPXfreeprob ((*p)->cplex_env, &((*p)->cplex_lp));
            }
            CPXcloseCPLEX (&((*p)->cplex_env));
        }
        COLOR_FREE (*p, COLORlp);
    }
}

int COLORlp_optimize (COLORlp *p)
{
    int rval = 0;
    int solstat;
    
    /** When restarting coloring with a set of stable sets (-r <filename>).  
        The Simplex Method has a very long running time for solving the first LP.
        On C4000.5 the dual simplex would take several hours.
        In such situations we presolve the first LP with the Barrier Method.

        The steepest edge norms are not initialized properly but only
        to 1. Therefore, Barrier should not be used when there are
        only a few rows.
     */
    if (p->noptcalls == 0) {
       int ncols = CPXgetnumcols (p->cplex_env, p->cplex_lp);
       int nrows = CPXgetnumrows (p->cplex_env, p->cplex_lp);
       if (ncols > 2 * nrows) {
          rval = CPXbaropt (p->cplex_env, p->cplex_lp);
          COLORcheck_rval (rval, "CPXbaropt failed");
       } 
    }
    rval = CPXdualopt (p->cplex_env, p->cplex_lp);
    COLORcheck_rval (rval, "CPXdualopt failed");
    

    solstat = CPXgetstat (p->cplex_env, p->cplex_lp);
    if (solstat == CPX_STAT_INFEASIBLE) {
        printf ("Infeasible LP\n");
        COLORlp_write (p, "joethelion.lp");
        rval = 2;  goto CLEANUP;
    } else if (solstat != CPX_STAT_OPTIMAL       &&
               solstat != CPX_STAT_OPTIMAL_INFEAS  ) {
        fprintf (stderr, "Cplex optimization status %d\n", solstat);
        if (solstat == CPX_STAT_ABORT_IT_LIM) {
            int itlim;
            rval = CPXgetintparam (p->cplex_env, CPX_PARAM_ITLIM, &itlim);
            if (!rval) {
                printf ("cplex iteration limit: %d\n", itlim);
                fflush (stdout);
            }
        }
        rval  = 1;  goto CLEANUP;
    }
    (p->noptcalls)++;

CLEANUP:
    return rval;
}

int COLORlp_objval (COLORlp *p, double *obj)
{
    int rval = 0;

    rval = CPXgetobjval (p->cplex_env, p->cplex_lp, obj);
    COLORcheck_rval (rval, "CPXgetobjval failed");

CLEANUP:
    return rval;
}

int COLORlp_change_objective(COLORlp *p, int start, int len, double* values)
{
    int i, rval = 0;
    int *indices = (int *) NULL;

    indices = COLOR_SAFE_MALLOC (len, int);
    COLORcheck_NULL (indices, "out of memory for indices");
    for (i = 0; i < len; i++) {
        indices[i] = start+i;
    }

    rval = CPXchgobj (p->cplex_env, p->cplex_lp, len, indices, values);
    COLORcheck_rval (rval, "CPXchgobj failed");

CLEANUP:
   COLOR_IFFREE (indices, int);
   return rval;
}

int COLORlp_addrow (COLORlp *p, int nzcount, int *cind, double *cval, 
       char sense, double rhs, char *name) 
{
    int rval = 0;
    char isense[1];
    char *iname[1];
    double irhs[1];
    int matbeg[1];


    switch (sense) {
    case COLORlp_EQUAL: 
        isense[0] = 'E'; break;
    case COLORlp_LESS_EQUAL:    
        isense[0] = 'L'; break;
    case COLORlp_GREATER_EQUAL:    
        isense[0] = 'G'; break;
    default: 
        fprintf (stderr, "unknown variable sense: %c\n", sense);
        rval = 1;  goto CLEANUP;
    }
 
    irhs[0] = rhs;
    iname[0] = name;
    matbeg[0] = 0;

    if (nzcount == 0) {
        rval = CPXnewrows (p->cplex_env, p->cplex_lp, 1, irhs,
                       isense, (double *) NULL, (char **) NULL);
        COLORcheck_rval (rval, "CPXnewrows failed");
    } else {
        rval = CPXaddrows (p->cplex_env, p->cplex_lp, 0, 1, nzcount, irhs,
                  isense, matbeg, cind, cval, (char **) NULL, (char **) NULL);
        COLORcheck_rval (rval, "CPXaddrows failed");
    }

CLEANUP:
    return rval;
}

int COLORlp_addcol (COLORlp *p, int nzcount, int *cind, double *cval, 
       double obj, double lb, double ub, char sense, char *name) 
{
    int rval = 0;
    int matbeg[1];
    double iobj[1], ilb[1], iub[1];
    char *iname[1];

    if (sense < 0) {
        printf ("bad sense, but cplex does not use this anyway\n");
        rval = 1; goto CLEANUP;
    }

    iobj[0] = obj;
    ilb[0] = lb;
    iub[0] = ub;
    iname[0] = name;
    matbeg[0] = 0;

    rval = CPXaddcols (p->cplex_env, p->cplex_lp, 1, nzcount, iobj, matbeg,
                       cind, cval, ilb, iub, (char **) NULL);
    COLORcheck_rval (rval, "CPXaddcols failed");

CLEANUP:
    return rval;
}

int COLORlp_deletecols (COLORlp *p, int first_cind, int last_cind)
{
    int rval = 0;

    rval = CPXdelcols (p->cplex_env, p->cplex_lp, first_cind, last_cind);
    COLORcheck_rval (rval, "CPXdelcols failed");

CLEANUP:
   return rval;
}


int COLORlp_pi (COLORlp *p, double *pi)
{
    int rval = 0;
    int nrows;

    nrows = CPXgetnumrows (p->cplex_env, p->cplex_lp);
    rval = CPXgetpi (p->cplex_env, p->cplex_lp, pi, 0, nrows - 1);
    COLORcheck_rval (rval, "CPXgetpi failed"); 

CLEANUP:
    return rval;
}

int COLORlp_x (COLORlp *p, double *x)
{
    int rval = 0;
    int ncols;

    ncols = CPXgetnumcols (p->cplex_env, p->cplex_lp);
    rval = CPXgetx (p->cplex_env, p->cplex_lp, x, 0, ncols - 1);
    COLORcheck_rval (rval, "CPXgetx failed");

CLEANUP:
    return rval;
}

int COLORlp_basis_cols (COLORlp *p, int *cstat)
{
   int rval = 0;
   int* rstat = (int*) NULL;

   rval =  CPXgetbase(p->cplex_env, p->cplex_lp, cstat, rstat);
   COLORcheck_rval (rval, "CPXgetbase failed");

 CLEANUP:
   return rval;
}

int COLORlp_set_all_coltypes (COLORlp *p, char sense)
{
   int rval = 0;

   if (!p) {
       printf ("COLORlp_set_all_coltypes called without an LP\n");
       rval = 1;  goto CLEANUP;
   }

   if (sense != COLORlp_CONTINUOUS) {
       printf ("Not set up to parse integer variables\n");
       rval = 1;  goto CLEANUP;
   }

CLEANUP:
   return rval;
}

int COLORlp_objective_sense (COLORlp *p, int sense)
{
    int rval = 0;
    char isense;

    if (sense == COLORlp_MIN) isense = CPX_MIN;
    else                      isense = CPX_MAX;

    CPXchgobjsen (p->cplex_env, p->cplex_lp, isense);

    return rval;
}

int COLORlp_setbound (COLORlp *p, int col, char lower_or_upper, double bnd)
{
    int rval = 0;
    int cindex[1];
    double bd[1];
    char lu[1];

    cindex[0] = col;
    lu[0] = lower_or_upper;
    bd[0] = bnd;

    rval = CPXchgbds (p->cplex_env, p->cplex_lp, 1, cindex, lu, bd);
    COLORcheck_rval (rval, "CPXchgbds failed");

CLEANUP:
    return rval;
}

int COLORlp_setnodelimit (COLORlp *p, int mip_node_limit)
{
    int rval = 0;

    if (!p || mip_node_limit < 1) {
        printf ("called with empty limit data, not set up in cplex\n");
        rval = 1;
    }
    return rval;
}

int COLORlp_write (COLORlp *p, const char *fname)
{
    int rval = 0;

    rval = CPXwriteprob (p->cplex_env, p->cplex_lp, fname, "RLP");
    COLORcheck_rval (rval, "CPXwriteprob failed");

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
    printf ("CPLEX error code: %d\n", c);
    fflush (stdout);
}

double COLORlp_int_tolerance ()
{
   return int_tolerance;
}
