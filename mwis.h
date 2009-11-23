#ifndef __MWIS_H
#define __MWIS_H

#include "color.h"


typedef struct _MWISenv      MWISenv;
typedef struct _MWISgrb_env  MWISgrb_env;
typedef struct _MWISls_env   MWISls_env;


int COLORstable_wrapper(MWISenv** env,
                        COLORset** newsets, int* nnewsets,
                        int ncount, int ecount, const int elist[], double nweights[]);

int COLORstable_LS(MWISls_env** env,
                   COLORset** newsets, int* nnewsets, int ncount, 
                   int ecount, const int elist[], double nweights[]);

int COLORstable_gurobi(MWISgrb_env** env,
                       COLORset** newsets, int* nnewsets,
                       int ncount, int ecount, const int elist[], double nweights[]);

int COLORstable_write_mps(const char*  filename,
                          int ncount, int ecount, const int elist[], double nweights[]);


int COLORstable_write_dimacs(const char*  filename,
                             int ncount, int ecount, const int elist[], double nweights[]);

int COLORstable_write_dimacs_clique(const char*  filename,
                                    int ncount, int ecount, const int elist[], double nweights[]);


int COLORstable_freeenv(MWISenv** env);
int COLORstable_free_ls_env(MWISls_env** env);
int COLORstable_free_grb_env(MWISgrb_env** env);

#endif
