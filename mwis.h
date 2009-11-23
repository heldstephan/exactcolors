#ifndef __MWIS_H
#define __MWIS_H

#include "graph.h"
#include "color.h"


typedef struct _MWISenv      MWISenv;
typedef struct _MWISgrb_env  MWISgrb_env;
typedef struct _MWISls_env   MWISls_env;

int COLORstable_initenv(MWISenv** env, const char* pname,
                        int write_mwis);

/* Fill newsets with independent sets of node weight larger than
   cutoff, or leave it empty, if no such independent set exists.
*/
int COLORstable_wrapper(MWISenv** env,
                        COLORset** newsets, int* nnewsets,
                        int ncount, int ecount, const int elist[], COLORNWT nweights[],
                        COLORNWT cutoff);

int COLORstable_LS(MWISls_env** env,
                   COLORset** newsets, int* nnewsets, int ncount,
                   int ecount, const int elist[], 
                   const COLORNWT nweights[],COLORNWT cutoff);

int COLORstable_gurobi(MWISgrb_env** env,
                       COLORset** newsets, int* nnewsets,
                       int ncount, int ecount, const int elist[], COLORNWT nweights[],
                       COLORNWT cutoff);

int COLORstable_write_mps(const char*  filename,
                          int ncount, int ecount, const int elist[], 
                          const COLORNWT nweights[],
                          COLORNWT cutoff);


/*
  Passes the node weights of type double from dbl_nweights to weights
  of type COLORNWT in nweights while scaling them by a the value
  scalef, which is determined on the fly and written to *scalef.
 */
int COLOR_double2COLORNWT(COLORNWT nweights[],
                          COLORNWT* scalef,
                          const double    dbl_nweights[],
                          int ncount);


/*
  Passes the node weights of type COLORNWT from nweights to weights
  of type double in dbl_nweights while dividing them by a the value
  divider.
*/
int COLOR_COLORNWT2double(double         dbl_nweights[],
                          const COLORNWT nweights[],
                          COLORNWT       divider,
                          int            ncount);

/** Write the stable set problem in dimacs format, adding the value of
    cutoff in the comment.
*/
int COLORstable_write_dimacs(const char*  filename,
                             int ncount, int ecount, const int elist[], 
                             const COLORNWT nweights[],
                             COLORNWT cutoff);


/** Write the stable set problem as a clique problem in the complement
    graph in dimacs format, adding the value of cutoff in the comment.
*/
int COLORstable_write_dimacs_clique(const char*  filename,
                                    int ncount, int ecount, const int elist[], 
                                    const COLORNWT nweights[],
                                    COLORNWT cutoff);


int COLORstable_freeenv(MWISenv** env);
int COLORstable_free_ls_env(MWISls_env** env);
int COLORstable_free_grb_env(MWISgrb_env** env);

#endif
