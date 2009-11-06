#ifndef __LP_H
#define __LP_H

#include <gurobi_c.h>

typedef struct COLORlp {
    GRBenv *env;
    GRBmodel *model;
} COLORlp;

int  COLORlp_init (COLORlp **p, const char *name);
void COLORlp_free (COLORlp **p);

#define COLORlp_CONTINUOUS GRB_CONTINUOUS
#define COLORlp_BINARY     GRB_BINARY
#define COLORlp_INTEGER    GRB_INTEGER

int COLORlp_optimize (COLORlp *p);
int COLORlp_objval (COLORlp *p, double *obj);
int COLORlp_pi (COLORlp *p, double *pi);


int COLORlp_addrow (COLORlp *p, int nzcount, int *cind, double *cval,                  char sense, double rhs, char *name);

int COLORlp_addcol (COLORlp *p, int nzcount, int *cind, double *cval,   
       double obj, double lb, double ub, char vartype, char *name);

int COLORlp_write (COLORlp *p, const char *fname);
void COLORlp_printerrorcode (int c);

#endif  /* __LP_H */
