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

#include<stdlib.h>
#include<stdio.h>
#include<limits.h>
#include<float.h>
#include<math.h>
#include<assert.h>
#include <string.h>

#include "graph.h"

#include "color.h"

#include "mwis_sewell/mwss_ext.h"
#include "mwis.h"

static int it = 0;
static const int max_ngreedy_fails = 1;
static const int max_ngreedy_switchoff_fails = 10;


/*
   Pending further empirical studies, this value should be set somewhere
   between 0.5 and 0.9.
*/
static const double high_density_threshold = 0.8;


/* int COLORstable_clique_enum(COLORset** newsets, int* nnewsets, int ncount, */
/*                             int ecount, const int elist[], COLORNWT nweights[], */
/*                             COLORNWT cutoff); */

struct _MWISenv{
   MWISgrb_env* grb_env;
   MWISls_env*  ls_env;
   const char*  pname;
   int          write_mwis;
   int          ngreedy_fails;
};


static
int COLORstable_max_weighted_node(COLORset** newsets, int* nnewsets, int ncount,
                                  COLORNWT nweights[],COLORNWT* objval)
{
   int rval = 0;
   int i;
   *objval = -1;

   (*newsets)  = COLOR_SAFE_MALLOC(1,COLORset);
   (*nnewsets) = 1;
   COLORcheck_NULL(*newsets,"Failed to allocate *newsets");

   (*newsets)[0].members = COLOR_SAFE_MALLOC(1,int);
   COLORcheck_NULL((*newsets)[0].members,"Failed to allocate (*newsets)[0].members");

   (*newsets)[0].count   = 1;

   for (i = 0; i < ncount;++i) {
      if (nweights[i] > *objval) {
         (*newsets)[0].members[0] = i;
         *objval = nweights[i];
      }
   }
 CLEANUP:
   return rval;
}

COLOR_MAYBE_UNUSED static
int COLORstable_clique_enum(COLORset** newsets, int* nnewsets, int ncount,
                            int ecount, const int elist[], COLORNWT nweights[],
                            COLORNWT cutoff, COLORNWT *objval)
{
   int            rval = 0;
   int            i;
   COLORadjgraph  G;
   COLORadjgraph  Gc;
   FILE*          file = (FILE*)  NULL;
   int*           elistc = (int*)  NULL;
   int            ecountc;
   COLORNWT*      oster_nweights = (COLORNWT*)  NULL;
   int            oster_cutoff = cutoff < COLOR_MAXINT ? cutoff + 1 : cutoff;

   rval = COLORadjgraph_build(&G,ncount,ecount,elist);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_build.");

   rval = COLORadjgraph_delete_unweighted(&G,&oster_nweights,nweights);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_delete_unweighted.");

   rval = COLORadjgraph_build_complement(&Gc, &G);
   COLORcheck_rval(rval,"Failed in COLORadjgraph_build_complement.");

   if (Gc.ecount) {
      COLORadjgraph_extract_edgelist(&ecountc,&elistc,&Gc);

      rval = COLORclique_ostergard(newsets,nnewsets,G.ncount,
                                   ecountc,elistc,oster_nweights,oster_cutoff,objval);
      COLORcheck_rval(rval,"Failed in COLORclique_enum.");
   } else {
      rval = COLORstable_max_weighted_node(newsets,nnewsets,G.ncount,
                                           oster_nweights,objval);
      COLORcheck_rval(rval,"Failed in COLORstable_max_weighted_node.");
   }
   if (COLORdbg_lvl() > 0) {
      printf("Best enumeration:   %13.10e ( %lld / %lld ).\n",
             COLORsafe_lower_dbl(*objval,cutoff),(long long ) *objval,(long long ) cutoff);
   }

   if (*objval > cutoff) {
      if (*nnewsets) {
         int oster_i = 0;
         int orig_i  = 0;
         while ( ! nweights[orig_i] ) {orig_i ++;}

         if (COLORdbg_lvl() > 0) {printf("NEW SET ");}
         
         for (i = 0; i < (*newsets)->count;++i) {
            while (oster_i < (*newsets)->members[i]) {
               oster_i++;
               orig_i++;
               while ( ! nweights[orig_i] ) {orig_i++;}
            }
            (*newsets)->members[i] = orig_i;
            if (COLORdbg_lvl() > 0) {printf(" %d",(*newsets)->members[i]);}
         }
         if (COLORdbg_lvl() > 0) {printf("\n");}

         COLORcheck_set((*newsets),ncount,ecount,elist);

/*          /\* Reverse order *\/ */
/*          for (i = 0; i < n/2 ;++i) { */
/*             int t; */
/*             COLOR_SWAP((*newsets)->members[i],(*newsets)->members[n-i-1],t); */
/*          } */
/*          printf("SWAPPED SET "); */
/*          for (i = 0; i < (*newsets)->count;++i) { */
/*             printf(" %d",(*newsets)->members[i]); */
/*          } */
/*          printf("\n"); */


      }
   } else {
      if (COLORdbg_lvl() > 0) {
         printf("BEST SET ");
         for (i = 0; i < (*newsets)->count;++i) {
            printf(" %d",(*newsets)->members[i]);
         }
         printf("\n");
      }
      COLORfree_sets(newsets,nnewsets);
   }


 CLEANUP:
   COLORadjgraph_free(&G);
   COLORadjgraph_free(&Gc);
   if (file) fclose(file);
   if (oster_nweights) free(oster_nweights);

   if (elistc) free(elistc);

   return rval;
}

static 
int COLORstable_sewell(COLORset** newsets, int* nnewsets, int ncount,
                       int ecount, const int elist[], COLORNWT nweights[],
                       COLORNWT cutoff,COLORNWT * objval)
{
   int            rval = 0;
   int            i;
   int            sewell_cutoff = cutoff < COLOR_MAXINT ? cutoff + 1 : cutoff;

   *objval   = .0;
   *newsets  = COLOR_SAFE_MALLOC(1,COLORset);
   *nnewsets = 1;

   (*newsets)->members = (int*) NULL;
   (*newsets)->count   = 0;

   
   rval = SEWELL_optimize( &((*newsets)->members),&((*newsets)->count),
                           ncount, ecount,elist,nweights,sewell_cutoff);
   COLORcheck_rval(rval,"Failed in SEWELL_optimize");

   qsort((*newsets)->members,(*newsets)->count,sizeof(int),
         COLORnode_comparator);

   for (i = 0; i < (*newsets)->count; ++i) {
      *objval += nweights[(*newsets)->members[i]];
   }

   if (COLORdbg_lvl() > 0) {
      printf("Best sewell:   %13.10e ( %lld / %lld ).\n",
             COLORsafe_lower_dbl(*objval,cutoff),(long long ) *objval,(long long ) cutoff);
   }

   if (*objval > cutoff) {
      if (*nnewsets) {
         if (COLORdbg_lvl() > 0) {
            printf("NEW SET ");
            for (i = 0; i < (*newsets)->count;++i) {
               printf(" %d",(*newsets)->members[i]);
            }
            printf("\n");
            
            COLORcheck_set((*newsets),ncount,ecount,elist);
         }
      }
   } else {
      if (COLORdbg_lvl() > 0) {
         printf("BEST SET ");
         for (i = 0; i < (*newsets)->count;++i) {
            printf(" %d",(*newsets)->members[i]);
         }
         printf("\n");
      }
      COLORfree_sets(newsets,nnewsets);
   }

 CLEANUP:
          
   return rval;
}

int COLORstable_initenv(MWISenv** env, const char* pname,
                        int write_mwis)
{
   int rval = 0;

   if (!*env) {
      *env = (MWISenv*) COLOR_SAFE_MALLOC (1,MWISenv);
      COLORcheck_NULL(*env,"Failed to allocate *env");
   }

   (*env)->grb_env = (MWISgrb_env*) NULL;
   (*env)->ls_env  = (MWISls_env*) NULL;
   (*env)->pname   = pname;
   (*env)->write_mwis   = write_mwis;

   (*env)->ngreedy_fails = 0;

 CLEANUP:
   return rval;
}

int COLORstable_wrapper(MWISenv** env,
                        COLORset** newsets, int* nnewsets, int ncount,
                        int ecount, const int elist[], COLORNWT nweights[],
                        COLORNWT cutoff)
{
   int rval = 0;
   double rtime;
   double density =  ((double) ecount) / ((double) (ncount * ( ncount - 1))) * 2.0;

   ++it;
   if (!*env) {
      const char* default_pname = "COLOR";
      int default_write_mwis    = 0;
      COLORstable_initenv(env,default_pname,default_write_mwis);
   }

   if (( *env)->ngreedy_fails ==  max_ngreedy_fails && COLORdbg_lvl() > 0) {
      printf("Greedy failed %d times in a row => not using greedy any more.\n",
             ( *env)->ngreedy_fails);
   }

   if ( (density < high_density_threshold) &&
         (( *env)->ngreedy_fails < max_ngreedy_switchoff_fails) ) {
      rtime = COLORcpu_time();
      rval = COLORstable_LS(&((*env)->ls_env),newsets, nnewsets,
                            ncount, ecount, elist,
                            nweights,cutoff);
      COLORcheck_rval(rval,"COLORstable_LS failed");
      rtime = COLORcpu_time() - rtime;
      if (COLORdbg_lvl() > 0) { printf("Greedy took %f seconds\n",rtime);}
   }

   if (*nnewsets) {
      ( *env)->ngreedy_fails = 0;
   } else {

      ++(( *env)->ngreedy_fails);

      if (!(*env)->ls_env) {
	 rval = COLORstable_init_LS(&((*env)->ls_env),
				    ncount,
				    ecount, elist,
                                       nweights,cutoff);
	 COLORcheck_rval(rval,"Failed in COLORstable_init_LS");
      }

      rval = COLORstable_round_down_weights((*env)->ls_env,
                                            nweights,cutoff);

      if( (*env)->write_mwis) {
         char filename [256];
         printf("Writing instances %s.*.%d.*\n",(*env)->pname,it);
         sprintf(filename,"%s.mwclq.%d.dimacs",(*env)->pname,it);
         COLORstable_write_dimacs_clique(filename,
                                         ncount, ecount, elist,
                                         nweights,cutoff);

         sprintf(filename,"%s.mwis.%d.dimacs",(*env)->pname,it);
         COLORstable_write_dimacs(filename,
                                  ncount, ecount, elist,
                                  nweights,cutoff);


         sprintf(filename,"%s.mwis.%d.lp",(*env)->pname,it);
#ifdef USE_GUROBI
         COLORstable_write_mps(filename,
                               ncount, ecount, elist,
                               nweights,cutoff);
#endif
      }

      if ( ncount <= SEWELL_node_limit() && density < high_density_threshold) {
         COLORNWT sewell_objval;

         rtime = COLORcpu_time();
         rval = COLORstable_sewell(newsets, nnewsets,
                                   ncount, ecount, elist,
                                   nweights,cutoff,&sewell_objval);
         COLORcheck_rval(rval,"COLORstable_LS failed");
         rtime = COLORcpu_time() - rtime;
         if (COLORdbg_lvl() > 0) { printf("Clique enumeration took %f seconds\n",rtime);}

      } else {
         COLORNWT oster_objval;

         rtime = COLORcpu_time();
         rval = COLORstable_clique_enum(newsets, nnewsets,
                                        ncount, ecount, elist,
                                        nweights,cutoff, &oster_objval);
         COLORcheck_rval(rval,"COLORstable_LS failed");
         rtime = COLORcpu_time() - rtime;
         if (COLORdbg_lvl() > 0) { printf("Clique enumeration took %f seconds\n",rtime);}
      }
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
#ifdef USE_GUROBI
      ls_rval  = COLORstable_free_grb_env(&(*env)->grb_env);
      if(grb_rval) printf("COLORstable_free_grb_env failed.");
#endif

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

   fprintf(file,"p edge %d %d\n",ncount, ecount);
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
   int            i,rval = 0;
   COLORadjgraph  G;
   COLORadjgraph  Gc;
   FILE*          file = (FILE*)  NULL;
   int*           elistc = (int*)  NULL;
   int            ecountc;
   int*           new_nweights = (int*)  NULL;

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

   fprintf(file,"p edge %d %d\n",Gc.ncount, ecountc);
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

   double dbl_scalef = COLORDBLmin(max_prec_dbl,max_mwiswt);

   dbl_scalef /= (double) ncount;

   for (i = 0; i < ncount;++i) {
      max_dbl_nweight =
         COLORDBLmax(max_dbl_nweight,dbl_nweights[i]);
   }
   dbl_scalef /= COLORDBLmax(1.0,max_dbl_nweight);
   dbl_scalef  = floor(dbl_scalef);
   *scalef  = (COLORNWT) dbl_scalef;

   for (i = 0; i < ncount;++i) {
      double weight = dbl_nweights[i] * dbl_scalef;
      assert(weight < (double) COLORNWT_MAX);
      nweights[i] = (COLORNWT) weight;
   }
   return 0;
}

#ifndef COMPILE_FOR_VALGRIND
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
#else
int COLOR_COLORNWT2double(double         dbl_nweights[],
                          const COLORNWT nweights[],
                          COLORNWT       divider,
                          int            ncount)
{
   int    i;
 
   double div_multiplier = (double) divider;
   div_multiplier = nextafter(div_multiplier,-DBL_MAX);
   div_multiplier = 1 / div_multiplier;
   div_multiplier = nextafter(div_multiplier, DBL_MAX);

   assert(divider > 0);

   for (i = 0; i < ncount;++i) {
      dbl_nweights[i] = ((double) nweights[i]) * div_multiplier  + DBL_EPSILON;
      dbl_nweights[i] = nextafter(dbl_nweights[i],DBL_MAX);
   }

   return 0;
}
#endif

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

   buf = (char*) COLOR_SAFE_MALLOC (bufsize, char);
   COLORcheck_NULL(buf,"Failed to allocate buf");

   setbuffer = (int*) COLOR_SAFE_MALLOC (ncount,int);
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
   *newsets = (COLORset*) COLOR_SAFE_MALLOC(*nnewsets,COLORset);
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
         int  unsorted = 0;
         token = strtok(p,delim);
         while( (token = strtok((char*) NULL,delim)) != (char*) NULL) {
            sscanf (token, "%d", &(setbuffer[setsize++]));
            if (setsize > 1 ){
               if (setbuffer[setsize-1] < setbuffer[setsize-2]) {
                  unsorted = 1;
               }
            }
         }
         if (unsorted) {
            qsort(setbuffer,setsize,sizeof(int),COLORnode_comparator);
         }
         (*newsets)[*nnewsets].count = setsize;
         (*newsets)[*nnewsets].members = (int*) COLOR_SAFE_MALLOC(setsize,int);
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

