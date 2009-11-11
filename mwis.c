#include<stdlib.h>
#include<stdio.h>
#include<limits.h>
#include<assert.h>
#include <sys/resource.h>
#include<time.h>

#include "graph.h"

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


int COLORstable_write_dimacs_clique(const char*  filename,
                                    int ncount, int ecount, const int elist[], double nweights[])
{
   int    rval = 0; 
   graph  G;
   graph  Gc;
   FILE*  file = (FILE*)  NULL;
   int*   elistc = (int*)  NULL;

   int*   old_nweights = (int*)  NULL;
   int*   new_nweights = (int*)  NULL;
      
   int    i,ecountc = 0;
   int    scalef;
   
   scalef = (double) (INT_MAX / (ncount + 1));
   
   old_nweights = (int*) malloc (ncount * sizeof(int));
   COLORcheck_NULL(old_nweights,"Failed to allocate old_nweights");
   
   for (i = 0; i < ncount; ++i) {
      old_nweights[i] = (int) ((double) scalef * nweights[i]);
   }

   rval = COLORadjgraph_build(&G,ncount,ecount,elist);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_build.");

   COLORadjgraph_delete_unweighted(&G,&new_nweights,old_nweights);
   ncount = G.ncount;
   


   rval = COLORadjgraph_build_complement(&Gc, &G);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_build_complement.");


   file = fopen(filename,"w");
   COLORcheck_NULL(file,"Failed to open clique file"); 

   COLORadjgraph_extract_edgelist(&ecountc,&elistc,&Gc);

   
   fprintf(file,"c Maximum weighted clique instances generated from exactcolor.\n");
   fprintf(file,"c scalef: %d.\n",scalef);
   fprintf(file,"A clique of value > scalef (%d) is a \n",scalef);

   fprintf(file,"p clq %d %d\n",Gc.ncount, ecountc);
   for (i = 0; i < ecountc;++i) {
      fprintf(file,"e %d %d\n",1+elistc[2*i],1+elistc[2*i+1]);
   }

   for (i = 0; i < ncount;++i) {
      assert(nweights[i] <= 1.0);
      fprintf(file,"n %d %d\n",i+1,new_nweights[i]);
   }


 CLEANUP:
   COLORadjgraph_free(&G);
   COLORadjgraph_free(&Gc);
   if (file) fclose(file);
   if (old_nweights) free(old_nweights);
   if (new_nweights) free(new_nweights);
   if (elistc) free(elistc);
   return rval;
}

