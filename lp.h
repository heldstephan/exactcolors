#ifndef __LP_H
#define __LP_H
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

#include <gurobi_c.h>

typedef struct COLORlp {
    GRBenv *env;
    GRBmodel *model;
} COLORlp;

typedef struct COLORlp_warmstart {
    int      rcount;
    int      ccount;
    int     *rstat;
    int     *cstat;
    double  *dnorm;
} COLORlp_warmstart;

int  COLORlp_init (COLORlp **p, const char *name);
void COLORlp_free (COLORlp **p);

#define COLORlp_CONTINUOUS GRB_CONTINUOUS
#define COLORlp_BINARY     GRB_BINARY
#define COLORlp_INTEGER    GRB_INTEGER

#define COLORlp_EQUAL         GRB_EQUAL
#define COLORlp_LESS_EQUAL    GRB_LESS_EQUAL
#define COLORlp_GREATER_EUQAL GRB_GREATER_EQUAL

int COLORlp_optimize (COLORlp *p);
int COLORlp_objval (COLORlp *p, double *obj);
int COLORlp_pi (COLORlp *p, double *pi);
int COLORlp_x (COLORlp *p, double *x);

int COLORlp_change_objective(COLORlp *p, int start, int len, double* values);

int COLORlp_addrow (COLORlp *p, int nzcount, int *cind, double *cval, 
                    char sense, double rhs, char *name);

int COLORlp_addcol (COLORlp *p, int nzcount, int *cind, double *cval,   
                    double obj, double lb, double ub, char vartype, char *name);

int COLORlp_deletecol (COLORlp *p, int cind);

int COLORlp_set_all_coltypes (COLORlp *p, char sense);

int COLORlp_objective_sense (COLORlp *p, int sense);
int COLORlp_setbound (COLORlp *p, int col, char lower_or_upper, double bnd);
int COLORlp_setnodelimit (COLORlp *p, int mip_node_limit);

int COLORlp_write (COLORlp *p, const char *fname);
void COLORlp_free_warmstart (COLORlp_warmstart **w);
void COLORlp_printerrorcode (int c);

double COLORlp_int_tolerance (void);

#endif  /* __LP_H */
