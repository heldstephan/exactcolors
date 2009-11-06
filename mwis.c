#include<stdio.h>
#include <sys/resource.h>
#include<time.h>
#include "mwis.h"

static double COLORwall_time (void)
{
    return (double) time (0);
}

static double COLORcpu_time (void)
{
    struct rusage ru;
    double t;

    getrusage (RUSAGE_SELF, &ru);

    t = ((double) ru.ru_utime.tv_sec) +
        ((double) ru.ru_utime.tv_usec) / 1000000.0;
    return t;
}


int COLORstable_wrapper(COLORset** newsets, int* nnewsets, int ncount, 
                        int ecount, const int elist[], double nweights[])
{
   int rval = 0;
   double rtime = COLORcpu_time();

   rval = COLORstable_LS(newsets, nnewsets, ncount, ecount, elist, nweights);
   COLORcheck_rval(rval,"COLORstable_LS failed");
   rtime = COLORcpu_time() - rtime;
   if (COLORdbg_lvl()) { printf("Greedy took %f seconds\n",rtime);}
   
/*    COLORfree_sets(newsets,nnewsets); */

   if (*nnewsets == 0) {
      rtime = COLORcpu_time();
      rval = COLORstable_gurobi(newsets, nnewsets, ncount, ecount, elist, nweights);
      COLORcheck_rval(rval,"COLORstable_LS failed");
      rtime = COLORcpu_time() - rtime;
      if (COLORdbg_lvl()) { printf("Gurobi took %f seconds\n",rtime);}
   }
 CLEANUP:
   return rval;
}

