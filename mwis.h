#ifndef __MWIS_H
#define __MWIS_H
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
                        COLORNWT cutoff, int greedy_only, int force_rounding);

int COLORstable_LS(MWISls_env** env,
                   COLORset** newsets, int* nnewsets, int ncount,
                   int ecount, const int elist[],
                   const COLORNWT nweights[],COLORNWT cutoff);

int COLORstable_init_LS(MWISls_env** env,
                        int ncount,
                        int ecount, const int elist[],
                        const COLORNWT nweights[], COLORNWT cutoff);

int COLORstable_gurobi(MWISgrb_env** env,
                       COLORset** newsets, int* nnewsets,
                       int ncount, int ecount, const int elist[], COLORNWT nweights[],
                       COLORNWT cutoff);

int COLORstable_write_mps(const char*  filename,
                          int ncount, int ecount, const int elist[],
                          const COLORNWT nweights[],
                          COLORNWT cutoff);

int COLORstable_read_stable_sets(COLORset** newsets, int* nnewsets,
                                 int ncount,
                                 const char* fname,
                                 const char* proplem_name);

int COLORstable_write_stable_sets(const COLORset* sets, int nsets,
                                  int ncount,
                                  const char* fname,
                                  const char* problem_name);

int COLORstable_round_down_weights(MWISls_env* env,
                                   COLORNWT nweights[],
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
