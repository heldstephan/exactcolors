#include<stdlib.h>
#include<stdio.h>
#include<limits.h>
#include<assert.h>
#include<math.h>
#include<float.h>

#include "graph.h"

#include "mwis.h"


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
      COLORcheck_NULL(*env,"Failed to allocate *env");

      (*env)->grb_env = (MWISgrb_env*) NULL;
      (*env)->ls_env = (MWISls_env*) NULL;
   }


   rval = COLORstable_LS(&((*env)->ls_env),newsets, nnewsets, ncount, ecount, elist, nweights);
   COLORcheck_rval(rval,"COLORstable_LS failed");
   rtime = COLORcpu_time() - rtime;
   if (COLORdbg_lvl() >= 0 ) { printf("Greedy took %f seconds\n",rtime);}
   
/*    COLORfree_sets(newsets,nnewsets); */

   if (*nnewsets == 0) {
      rtime = COLORcpu_time();
      rval = COLORstable_gurobi(&((*env)->grb_env),newsets, nnewsets, ncount, ecount, elist, nweights);
      COLORcheck_rval(rval,"COLORstable_LS failed");
      rtime = COLORcpu_time() - rtime;
      if (COLORdbg_lvl() >= 0) { printf("Gurobi took %f seconds\n",rtime);}
   }

 CLEANUP:
   if (rval) {
      COLORstable_freeenv(env);
   }
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

int COLORstable_write_dimacs(const char*  filename,
                             int ncount, 
                             int ecount, 
                             const int elist[], 
                             double nweights[])
{
   int    rval = 0; 
   FILE*  file = (FILE*)  NULL;
      
   int    i;
   
   int       lscalef  = (INT_MAX) / ncount;
   double     scalef  = (double) lscalef;



   file = fopen(filename,"w");
   COLORcheck_NULL(file,"Failed to open clique file"); 
   
   fprintf(file,"c Maximum weighted independent set instances generated from exactcolors.\n");
   fprintf(file,"c scalef: %d (== (INT_MAX) / ncount).\n",lscalef);
   fprintf(file,"c An independent set of value > %d defines an improving stable sets for coloring.\n",lscalef);

   fprintf(file,"p graph %d %d\n",ncount, ecount);
   for (i = 0; i < ecount;++i) {
      fprintf(file,"e %d %d\n",1+elist[2*i],1+elist[2*i+1]);
   }

   for (i = 0; i < ncount;++i) {
      assert(nweights[i] <= 1.0);
      fprintf(file,"n %d %llu\n",i+1,(long long) (scalef * nweights[i]));
   }

 CLEANUP:
   if (file) fclose(file);
   return rval;
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
   double scalef  = exp2(DBL_MANT_DIG-1) / ncount;
   
   scalef = (double) (LLONG_MAX / (ncount + 1));
   
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

   
   fprintf(file,"c Maximum weighted clique instances generated from exactcolors.\n");
   fprintf(file,"c scalef: %u.\n",(unsigned int)scalef);
   fprintf(file,"c A clique of value > scalef (%u) defines an improving stable sets for coloring.\n",(unsigned int) scalef);

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

