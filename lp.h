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

typedef struct COLORlp COLORlp;

typedef struct COLORlp_warmstart {
    int      rcount;
    int      ccount;
    int     *rstat;
    int     *cstat;
    double  *dnorm;
} COLORlp_warmstart;

int  COLORlp_init (COLORlp **p, const char *name);
void COLORlp_free (COLORlp **p);

int  COLORlp_init_env(void);
void COLORlp_free_env(void);

#define COLORlp_CONTINUOUS 0
#define COLORlp_BINARY     1
#define COLORlp_INTEGER    2

#define COLORlp_EQUAL         'E'
#define COLORlp_LESS_EQUAL    'L'
#define COLORlp_GREATER_EQUAL 'G'

#define COLORlp_LOWER      0
#define COLORlp_BASIC      1
#define COLORlp_UPPER      2
#define COLORlp_FREE       3

#define COLORlp_MIN  1
#define COLORlp_MAX -1

int COLORlp_optimize (COLORlp *p);
int COLORlp_objval (COLORlp *p, double *obj);
int COLORlp_pi (COLORlp *p, double *pi);
int COLORlp_x (COLORlp *p, double *x);

int COLORlp_basis_cols (COLORlp *p, int *cstat);

int COLORlp_change_objective(COLORlp *p, int start, int len, double* values);

int COLORlp_addrow (COLORlp *p, int nzcount, int *cind, double *cval,
                    char sense, double rhs, char *name);

int COLORlp_addcol (COLORlp *p, int nzcount, int *rind, double *rval,
                    double obj, double lb, double ub, char vartype, char *name);

int COLORlp_deletecols (COLORlp *p, int first_cind, int last_cind);

int COLORlp_set_all_coltypes (COLORlp *p, char sense);

int COLORlp_objective_sense (COLORlp *p, int sense);
int COLORlp_setbound (COLORlp *p, int col, char lower_or_upper, double bnd);
int COLORlp_setnodelimit (COLORlp *p, int mip_node_limit);
int COLORlp_set_cutoff (COLORlp *p, double cutoff);


int COLORlp_write (COLORlp *p, const char *fname);
void COLORlp_free_warmstart (COLORlp_warmstart **w);
void COLORlp_printerrorcode (int c);

double COLORlp_int_tolerance (void);

#endif  /* __LP_H */
