#ifndef __MWIS_H
#define __MWIS_H

#include "color.h"

int COLORstable_wrapper(COLORset** newsets, int* nnewsets,
                        int ncount, int ecount, const int elist[], double nweights[]);

int COLORstable_LS(COLORset** newsets, int* nnewsets, int ncount, 
                   int ecount, const int elist[], double nweights[]);

int COLORstable_gurobi(COLORset** newsets, int* nnewsets,
        int ncount, int ecount, const int elist[], double nweights[]);


#endif
