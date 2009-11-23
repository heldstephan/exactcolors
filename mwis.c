#include<stdlib.h>
#include<stdio.h>
#include<limits.h>
#include<float.h>
#include<math.h>
#include<assert.h>
#include <string.h>

#include "graph.h"

#include "mwis.h"


struct _MWISenv{
   MWISgrb_env* grb_env;
   MWISls_env*  ls_env;
   const char*  pname;
   int          write_mwis;
};

static int it = 0;

int COLORstable_initenv(MWISenv** env, const char* pname,
                        int write_mwis)
{
   int rval = 0;

   if (!*env) {
      *env = (MWISenv*) malloc(sizeof(MWISenv));
      COLORcheck_NULL(*env,"Failed to allocate *env");
   }
   
   (*env)->grb_env = (MWISgrb_env*) NULL;
   (*env)->ls_env  = (MWISls_env*) NULL;
   (*env)->pname   = pname;
   (*env)->write_mwis   = write_mwis;

 CLEANUP:
   return rval;
}

int COLORstable_wrapper(MWISenv** env,
                        COLORset** newsets, int* nnewsets, int ncount,
                        int ecount, const int elist[], COLORNWT nweights[],
                        COLORNWT cutoff)
{
   int rval = 0;
   double rtime = COLORcpu_time();
   
   ++it;
   if (!*env) {
      const char* default_pname = "COLOR";
      int default_write_mwis    = 0;
      COLORstable_initenv(env,default_pname,default_write_mwis);
   }


   rval = COLORstable_LS(&((*env)->ls_env),newsets, nnewsets, 
                         ncount, ecount, elist, 
                         nweights,cutoff);
   COLORcheck_rval(rval,"COLORstable_LS failed");
   rtime = COLORcpu_time() - rtime;
   if (COLORdbg_lvl() >= 0 ) { printf("Greedy took %f seconds\n",rtime);}

   /* Uncomment to enforce gurobi.*/
   /*    COLORfree_sets(newsets,nnewsets); */

   if (*nnewsets == 0) {
      if( (*env)->write_mwis) {
         char filename [256];
         sprintf(filename,"%s.mwclq.%d.dimacs",(*env)->pname,it);
         COLORstable_write_dimacs_clique(filename,
                                         ncount, ecount, elist, 
                                         nweights,cutoff);
         
         sprintf(filename,"%s.mwis.%d.dimacs",(*env)->pname,it);
         COLORstable_write_dimacs(filename,
                                  ncount, ecount, elist, 
                                  nweights,cutoff);
         
         
         sprintf(filename,"%s.mwis.%d.lp",(*env)->pname,it);
         COLORstable_write_mps(filename,
                               ncount, ecount, elist, 
                               nweights,cutoff);
      }

      rtime = COLORcpu_time();
      rval = COLORstable_gurobi(&((*env)->grb_env),newsets, nnewsets, 
                                ncount, ecount, elist, 
                                nweights,cutoff);
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
                             const COLORNWT nweights[],
                             COLORNWT cutoff)
{
   int    rval = 0;
   FILE*  file = (FILE*)  NULL;
   int    i;


   file = fopen(filename,"w");
   COLORcheck_NULL(file,"Failed to open mwis file");
   
   fprintf(file,"c Maximum weighted independent set instances generated from exactcolors.\n");
   fprintf(file,"c scalef: %d .\n",cutoff);
   fprintf(file,"c An independent set of value > %d defines an improving stable sets for coloring.\n",cutoff);

   fprintf(file,"p graph %d %d\n",ncount, ecount);
   for (i = 0; i < ecount;++i) {
      fprintf(file,"e %d %d\n",1+elist[2*i],1+elist[2*i+1]);
   }
   
   for (i = 0; i < ncount;++i) {
      assert(nweights[i] <= cutoff);
      fprintf(file,"n %d %llu\n",i+1,(long long) nweights[i]);
   }

 CLEANUP:
   if (file) fclose(file);
   return rval;
}


int COLORstable_write_dimacs_clique(const char*  filename,
                                    int ncount, int ecount, const int elist[], 
                                    const COLORNWT nweights[], COLORNWT cutoff)
{
   int    i,rval = 0;
   graph  G;
   graph  Gc;
   FILE*  file = (FILE*)  NULL;
   int*   elistc = (int*)  NULL;
   int    ecountc;
   int*   new_nweights = (int*)  NULL;

   rval = COLORadjgraph_build(&G,ncount,ecount,elist);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_build.");
   
   COLORadjgraph_delete_unweighted(&G,&new_nweights,nweights);
   ncount = G.ncount;


   rval = COLORadjgraph_build_complement(&Gc, &G);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_build_complement.");
   

   file = fopen(filename,"w");
   COLORcheck_NULL(file,"Failed to open clique file");
   
   COLORadjgraph_extract_edgelist(&ecountc,&elistc,&Gc);


   fprintf(file,"c Maximum weighted clique instances generated from "
           "exactcolors.\n");
   fprintf(file,"c cutoff: %lld.\n",(long long) cutoff);
   fprintf(file,"c A clique of value > scalef (%lld) defines an "
           "improving stable sets for coloring.\n",
           (long long) cutoff);
   
   fprintf(file,"p clq %d %d\n",Gc.ncount, ecountc);
   for (i = 0; i < ecountc;++i) {
      fprintf(file,"e %d %d\n",1+elistc[2*i],1+elistc[2*i+1]);
   }

   for (i = 0; i < ncount;++i) {
      assert(nweights[i] <= cutoff);
      fprintf(file,"n %d %d\n",i+1,new_nweights[i]);
   }


 CLEANUP:
   COLORadjgraph_free(&G);
   COLORadjgraph_free(&Gc);
   if (file) fclose(file);
   if (new_nweights) free(new_nweights);
   if (elistc) free(elistc);
   return rval;
}

int COLOR_double2COLORNWT(COLORNWT nweights[],
                          COLORNWT* scalef,
                          const double    dbl_nweights[],
                          int ncount)
{
   int    i;
   double max_dbl_nweight= -DBL_MAX;
   double max_prec_dbl = exp2(DBL_MANT_DIG-1);
   static const double max_mwiswt   = (double) COLORNWT_MAX;

   double dbl_scalef = dblmin(max_prec_dbl,max_mwiswt);

   dbl_scalef /= (double) ncount;

   for (i = 0; i < ncount;++i) {
      max_dbl_nweight =
         dblmax(max_dbl_nweight,dbl_nweights[i]);
   }
   dbl_scalef /= dblmax(1.0,max_dbl_nweight);
   dbl_scalef  = floor(dbl_scalef);
   *scalef  = (COLORNWT) dbl_scalef;

   for (i = 0; i < ncount;++i) {
      double weight = dbl_nweights[i] * dbl_scalef;
      assert(weight < (double) COLORNWT_MAX);
      nweights[i] = (COLORNWT) weight;
   }
   return 0;
}

int COLOR_COLORNWT2double(double         dbl_nweights[],
                          const COLORNWT nweights[],
                          COLORNWT       divider,
                          int            ncount)
{
   int    i;
   int current_rounding = fegetround();

   fesetround(FE_DOWNWARD);   
   double div_multiplier = (double) divider;

   fesetround(FE_UPWARD);  
   div_multiplier = 1 / div_multiplier;

   assert(divider > 0);

   for (i = 0; i < ncount;++i) {
      dbl_nweights[i] = ((double) nweights[i]) * div_multiplier + DBL_EPSILON;
   }

   fesetround(current_rounding);
   
   return 0;
}

int COLORstable_read_stable_sets(COLORset** newsets, int* nnewsets,
                                 int ncount,
                                 const char* fname,
                                 const char* problem_name)
{
   int rval = 0;

   FILE* ifile = (FILE*) NULL;    
   char* buf = (char*) NULL;
   char* p;
   int  bufsize = 2 * ncount * (2 + (int) ceil(log((double)ncount + 10)));
   int* setbuffer = (int*) NULL;
   const char* delim = " \t\n";
   char* token = (char* ) NULL;

   buf = (char*) malloc(bufsize);
   COLORcheck_NULL(buf,"Failed to allocate buf");

   setbuffer = (int*) malloc(ncount * sizeof(int));
   COLORcheck_NULL(setbuffer,"Failed to allocate setbuffer");
   
   ifile = fopen(fname,"r");
   if(!ifile) {
      fprintf(stderr, "Failed to open %s for reading",fname);
      rval = 1; goto CLEANUP;
   }

   *nnewsets = 0;
   while (fgets (buf, bufsize, ifile) != (char *) NULL) {
      p = buf;
      if (p[0] == 's') {
         (*nnewsets)++;
      }
   }
   *newsets = (COLORset*) malloc(*nnewsets * sizeof(COLORset));
   COLORcheck_NULL(newsets,"Failed to allocate *newsets");
   while( *nnewsets > 0) {
      (*nnewsets)--;
      (*newsets)[*nnewsets].members = (int*) NULL;
   }

   

   rewind(ifile);
   while (fgets (buf, bufsize, ifile) != (char *) NULL) {
      p = buf;
      if (p[0] == 'c') {
         printf ("Comment: %s", p+1);
      } else if (p[0] == 'p') {
         token = strtok(p,delim); /* get 'p' */ 
         
         token = strtok((char*) NULL,delim); /* get problem name */
         if(strcmp(token,problem_name)) {
            fprintf(stderr,"Stable set file problem %s does not match instance problem %s",
                    token, problem_name);
            rval = 1; goto CLEANUP;
         }
         token = strtok(NULL,delim);

      } else if( p[0] == 's') {
         int  setsize = 0;
         token = strtok(p,delim);
         while( (token = strtok((char*) NULL,delim)) != (char*) NULL) {
            sscanf (token, "%d", &(setbuffer[setsize++]));
         }
         (*newsets)[*nnewsets].count = setsize;
         (*newsets)[*nnewsets].members = (int*) malloc(setsize * sizeof(int));
         COLORcheck_NULL((*newsets)[*nnewsets].members,
                         "Failed to allocate (*newsets)[*nnewsets].members");
         memcpy((*newsets)[*nnewsets].members, setbuffer,setsize * sizeof(int));
         (*nnewsets)++;
      }
   }

   if (! (*nnewsets)) {
      fprintf(stderr,"Failed to read any sets");
      rval = 1; goto CLEANUP;
   }
   
 CLEANUP:
   if (ifile) fclose(ifile);
   return rval; 
}

int COLORstable_write_stable_sets(const COLORset* sets, int nsets,
                                  int   ncount,
                                  const char* fname,
                                  const char* problem_name)
{
   int rval = 0;
   int i;
   FILE* ofile = (FILE*) NULL;
   
   ofile = fopen(fname,"w");
   if(!ofile) {
      fprintf(stderr, "Failed to open %s for writing",fname);
      rval = 1; goto CLEANUP;
   }
   

   fprintf(ofile,"c Stable sets from exactcolors on instance %s:\n",problem_name);
   fprintf(ofile,"c Syntax:\n");
   fprintf(ofile,"c \t p <problem name> <node count>\n");
   fprintf(ofile,"c \t s <list of intergers ( = nodes of a stable set)>\n");

   fprintf(ofile,"p %s %d\n",problem_name,ncount);
   for(i = 0; i < nsets; ++i) {
      int j;
      int  count   = sets[i].count;
      int* members = sets[i].members;
      fprintf(ofile,"s");
      for(j = 0; j < count; ++j) {
         fprintf(ofile," %d",members[j]);
      }
      fprintf(ofile,"\n");
   }
   
 CLEANUP:
   if (ofile) fclose(ofile);
   return rval; 
}

