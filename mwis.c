#include<stdio.h>
#include "mwis.h"

int COLORstable_wrapper(COLORset** newsets, int* nnewsets, int ncount, 
                        int ecount, const int elist[], double nweights[])
{
   int rval = 0;
   rval = COLORstable_LS(newsets, nnewsets, ncount, ecount, elist, nweights);
   COLORcheck_rval(rval,"COLORstable_LS failed");
   
   if (*nnewsets == 0) {
      rval = COLORstable_gurobi(newsets, nnewsets, ncount, ecount, elist, nweights);
      COLORcheck_rval(rval,"COLORstable_LS failed");
   }
 CLEANUP:
   return rval;
}

