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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lp.h"
#include "color.h"
#include <qsopt.h>

struct COLORlp {
    QSprob p;
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

    (*p)->p = QScreate_prob (name, QS_MIN);
    if ((*p)->p == (QSprob) NULL) {
        fprintf (stderr, "QScreate_prob failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSset_param ((*p)->p, QS_PARAM_DUAL_PRICING, QS_PRICE_DSTEEP);
    if (rval) {
        fprintf (stderr, "QSset_param failed\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:
    return rval;
}

void COLORlp_free (COLORlp **p)
{
    if (*p) {
        QSfree_prob ((*p)->p);
        COLOR_FREE (*p, COLORlp);
    }
}

int COLORlp_optimize (COLORlp *p)
{
    int rval = 0;
    int status;

    rval = QSopt_dual (p->p, &status);
    if (rval) {
        fprintf (stderr, "QSopt_dual failed\n"); goto CLEANUP;
    }

    if (status == QS_LP_ITER_LIMIT) {
        printf ("Dual LP Solver reached iteration limit\n"); fflush (stdout);
        rval = QSwrite_prob (p->p, "iter.lp", "LP");
        if (rval) {
            fprintf (stderr, "QSwrite_prob failed\n");
            goto CLEANUP;
        }
        printf ("Saved LP as iter.lp\n"); fflush (stdout);
        rval = QSwrite_basis (p->p, (QSbas) NULL, "iter.bas");
        if (rval) {
            fprintf (stderr, "QSwrite_basis failed\n"); goto CLEANUP;
        }
        printf ("Saved LP as iter.bas\n"); fflush (stdout);
    } else if (status == QS_LP_TIME_LIMIT) {
        printf ("Dual LP Solver reached time limit\n"); fflush (stdout);
    }

    if (status == QS_LP_INFEASIBLE) {
        rval = 2; goto CLEANUP;
    } else if (status != QS_LP_OPTIMAL) {
        fprintf (stderr, "no optimal LP-solution exists: %d\n", status);
        rval = 1; goto CLEANUP;
    }

CLEANUP:
    return rval;
}

int COLORlp_objval (COLORlp *p, double *obj)
{
    int rval = 0;

    rval = QSget_objval (p->p, obj);
    if (rval) {
        fprintf (stderr, "QSget_objval failed");
        goto CLEANUP;
    }


CLEANUP:
    return rval;

}

int COLORlp_change_objective(COLORlp *p, int start, int len, double* values)
{
    int i, rval = 0;

    for (i = 0; i < len; i++) {
        rval = QSchange_objcoef (p->p, start+i, values[i]);
        COLORcheck_rval (rval, "Failed in QSchange_objcoef");
    }

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
        isense = 'E'; break;
    case COLORlp_LESS_EQUAL:
        isense = 'L'; break;
    case COLORlp_GREATER_EQUAL:
        isense = 'G'; break;
    default:
        fprintf (stderr, "unknown variable sense: %c\n", sense);
        rval = 1;  goto CLEANUP;
    }

    rval = QSadd_row (p->p, nzcount, cind, cval, rhs, isense, name);
    COLORcheck_rval (rval, "QSadd_row failed");

CLEANUP:
    return rval;
}

int COLORlp_addcol (COLORlp *p, int nzcount, int *cind, double *cval,
       double obj, double lb, double ub, char sense, char *name)
{
    int rval = 0;

    if (sense < 0) {
        printf ("bad sense, but qsopt does not use this anyway\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSadd_col (p->p, nzcount, cind, cval, obj, lb, ub, name);
    COLORcheck_rval (rval, "QSadd_col failed");

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

   rval = QSdelete_cols (p->p, numdel,dellist);
   COLORcheck_rval (rval, "QSdelete_col failed");

CLEANUP:
   if (dellist) free(dellist);
   return rval;
}


int COLORlp_pi (COLORlp *p, double *pi)
{
    int rval = 0;
    int status;

    rval = QSget_status (p->p, &status);
    if (rval) {
        fprintf (stderr, "QSget_status failed"); goto CLEANUP;
    }

    if (status == QS_LP_INFEASIBLE) {
        rval = QSget_infeas_array (p->p, pi);
        COLORcheck_rval (rval, "QSget_infeas_array failed");
    } else {
        rval = QSget_pi_array (p->p, pi);
        COLORcheck_rval (rval, "QSget_pi_array failed");
    }

CLEANUP:
    return rval;
}

int COLORlp_x (COLORlp *p, double *x)
{
    int rval = 0;

    rval = QSget_x_array (p->p, x);
    COLORcheck_rval (rval, "QSget_x_array failed");

 CLEANUP:
    return rval;
}

int COLORlp_basis_cols (COLORlp *p, int *int_cstat)
{
   int rval = 0;
   char* rstat = (char*) NULL;
   char* cstat = (char*) NULL;
   int   ncols = QSget_colcount (p->p);
   int   nrows = QSget_rowcount (p->p);
   int   i;
   cstat = COLOR_SAFE_MALLOC(ncols, char);
   COLORcheck_NULL(cstat,"Failed to allocate cstat");

   rstat = COLOR_SAFE_MALLOC(nrows, char);
   COLORcheck_NULL(rstat,"Failed to allocate rstat");


   rval = QSget_basis_array (p->p, cstat,rstat);
   COLORcheck_rval (rval, "QSget_basis_array failed");

   for (i = 0; i < ncols; ++i) {
      switch (cstat[i]) {
      case QS_COL_BSTAT_LOWER:
         int_cstat[i] = COLORlp_LOWER;
         break;
      case QS_COL_BSTAT_UPPER:
         int_cstat[i] = COLORlp_UPPER;
         break;
      case QS_COL_BSTAT_FREE:
         int_cstat[i] = COLORlp_FREE;
         break;
      case QS_COL_BSTAT_BASIC:
         int_cstat[i] = COLORlp_BASIC;
         break;
      default:
         rval = 1;
         COLORcheck_rval(rval,"ERROR: Received unknown cstat");
      }
   }
 CLEANUP:
   if(cstat) free(cstat);
   if(rstat) free(rstat);

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
       printf ("QSopt does not handled integer variables\n");
       rval = 1;  goto CLEANUP;
   }

CLEANUP:
   return rval;
}

int COLORlp_objective_sense (COLORlp *p, int sense)
{
    int rval = 0;
     char isense;

    if (sense == COLORlp_MIN) isense = QS_MIN;
    else                      isense = QS_MAX;

    rval = QSchange_objsense (p->p, isense);
    COLORcheck_rval (rval, "QSchange_objsense failed");

CLEANUP:
    return rval;
}

int COLORlp_setbound (COLORlp *p, int col, char lower_or_upper, double bnd)
{
    int rval = 0;

    rval = QSchange_bound (p->p, col, lower_or_upper, bnd);
    COLORcheck_rval (rval, "QSchange_bounds failed");

CLEANUP:
    return rval;
}

int COLORlp_setnodelimit (COLORlp *p, int mip_node_limit)
{
    int rval = 0;

    if (!p || mip_node_limit < 1) {
        printf ("called with empty limit data, no meaning in QSopt\n");
        rval = 1;
    }
    return rval;
}

int COLORlp_write (COLORlp *p, const char *fname)
{
    int rval = 0;

    rval = QSwrite_prob (p->p, fname, "LP");
    COLORcheck_rval (rval, "QSwrite_prob failed");

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
    printf ("QSopt error code: %d\n", c);
    fflush (stdout);
}

double COLORlp_int_tolerance ()
{
   return int_tolerance;
}

int COLORlp_set_cutoff (COLORlp *p, double cutoff)
{
   int rval = 1;

   (void) p;
   (void) cutoff;

   COLORcheck_rval(rval,"COLORlp_set_cutoff not yet implemented.");

CLEANUP:
   return rval;
}
