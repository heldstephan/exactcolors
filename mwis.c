#include<stdlib.h>
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


struct _MWISenv{
   MWISgrb_env* grb_env;
   MWISls_env*  ls_env;
};

int COLORstable_wrapper(MWISenv** env,
                        COLORset** newsets, int* nnewsets, int ncount, 
                        int ecount, const int elist[], double nweights[])
{
   int rval = 0;
   double rtime = COLORcpu_time();

   if (!*env) {
      *env = (MWISenv*) malloc(sizeof(MWISenv));
      (*env)->grb_env = (MWISgrb_env*) NULL;
      (*env)->ls_env = (MWISls_env*) NULL;
   }


/*    rval = COLORstable_LS(&((*env)->ls_env),newsets, nnewsets, ncount, ecount, elist, nweights); */
/*    COLORcheck_rval(rval,"COLORstable_LS failed"); */
/*    rtime = COLORcpu_time() - rtime; */
/*    if (COLORdbg_lvl()) { printf("Greedy took %f seconds\n",rtime);} */
   
/*    COLORfree_sets(newsets,nnewsets); */

   if (*nnewsets == 0) {
      rtime = COLORcpu_time();
      rval = COLORstable_gurobi(&((*env)->grb_env),newsets, nnewsets, ncount, ecount, elist, nweights);
      COLORcheck_rval(rval,"COLORstable_LS failed");
      rtime = COLORcpu_time() - rtime;
      if (COLORdbg_lvl()) { printf("Gurobi took %f seconds\n",rtime);}
   }
 CLEANUP:
   return rval;
}

int COLORstable_freeenv(MWISenv** env)
{
   int grb_rval = 0;
   int ls_rval = 0;
   if (*env) {
      grb_rval = COLORstable_free_ls_env(&(*env)->ls_env);
      if(ls_rval) printf("COLORstable_free_ls_env failed.");
      ls_rval  = COLORstable_free_grb_env(&(*env)->grb_env);
      if(grb_rval) printf("COLORstable_free_grb_env failed.");
      
      free(*env);
      *env = (MWISenv*) NULL;
   }
   return grb_rval + ls_rval;
}
