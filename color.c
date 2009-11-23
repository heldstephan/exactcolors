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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <getopt.h>
#include <string.h>
#include <fenv.h>
#include <sys/resource.h>
#include <time.h>

#include "lp.h"
#include "graph.h"
#include "color.h"

#include "mwis.h"
#include "plotting.h"

static char *edgefile = (char *) NULL;
static char *outfile = (char *) NULL;
static char *cclasses_outfile = (char *) NULL;
static char *cclasses_infile = (char *) NULL;


static int debug = 0;
static int write_mwis = 0;

int main (int ac, char **av);
static int parseargs (int ac, char **av);
static void usage (char *f);
static int get_problem_name(char* pname,const char* efname);

int COLORdbg_lvl() {
   return debug;
}

static void print_objective(COLORNWT*       lower_bound,
                            const COLORNWT* mwis_pi,
                            COLORNWT        mwis_pi_scalef,
                            int             ncount)
{
   int i;
   *lower_bound = .0;
   for (i = 0; i < ncount;++i) {
      *lower_bound += (double) mwis_pi[i];
   }

   printf("Current primal LP objective: %f.\n",
          COLORsafe_lower_dbl(*lower_bound,mwis_pi_scalef));
}

static void make_pi_feasible(double* pi,COLORset* gcolors,int gcount)
{
   int c;
   int current_rounding = fegetround();
   fesetround(FE_UPWARD);

   for (c = 0; c < gcount;++c) {
      int i;
      double colsum = .0;
      double newcolsum = .0;
      for (i = 0; i < gcolors[c].count;++i) {
         if (signbit(pi[gcolors[c].members[i]])) {
            pi[gcolors[c].members[i]] = 0.0;
         }
         colsum += pi[gcolors[c].members[i]];
      }
      if (colsum > 1.0) {
         fesetround(FE_DOWNWARD);
         
         for (i = 0; i < gcolors[c].count;++i) {
            pi[gcolors[c].members[i]] /= colsum;
            newcolsum += pi[gcolors[c].members[i]];
         }
         if (COLORdbg_lvl()> 1) {printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n",c,colsum,newcolsum);}
         fesetround(FE_UPWARD);
      }
      fesetround(current_rounding);
   }
}

   static int concat_newsets(COLORset** gcolors, int* gcount, int* gallocated,
                          COLORset** newsets, int* nnewsets)
   {
   int rval = 0;
   COLORset* tmpsets = (COLORset*) NULL;
   
   if (*nnewsets == 0) return rval;
   
   if (*gcount + *nnewsets > *gallocated) {
      /* Double size */
      tmpsets = malloc(2 * (*gallocated) * sizeof(COLORset));
      COLORcheck_NULL (tmpsets, "out of memory for tmpsets");
      memcpy(tmpsets,*gcolors, (*gcount) * sizeof(COLORset));
      free(*gcolors);
      *gallocated *= 2;
      *gcolors = tmpsets;
      tmpsets = NULL;
   }
   memcpy(*gcolors + *gcount, *newsets, (*nnewsets) * sizeof(COLORset));
   *gcount +=  *nnewsets;
 CLEANUP:
   if (rval) {
      if (*gcolors) free(*gcolors);
   }
   if (*newsets) free(*newsets);
   *newsets = (COLORset*) NULL;
   *nnewsets = 0;
   return rval;
}

int main (int ac, char **av)
{
    char pname[256] = "";
    int rval = 0;
    int iterations = 0;
    int maxiterations = 1000000;
    int ncount = 0, ecount = 0, gcount, gallocated, i, j;
    int *elist = (int *) NULL;
    double *coef = (double *) NULL;
    double *pi = (double *) NULL;
    COLORNWT  mwis_pi_scalef;
    COLORNWT *mwis_pi = (COLORNWT *) NULL;
    COLORNWT  lower_bound = 0;
    COLORNWT  lower_scaled_bound = 0;
    COLORset *gcolors = (COLORset *) NULL;
    COLORset *newsets = (COLORset*) NULL;   
    COLORlp *lp       = (COLORlp *) NULL;
    int       nnewsets = 0;
    int break_while_loop = 1;
    MWISenv* mwis_env  = (MWISenv*) NULL;
    double tot_rtime;
   
    
    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    get_problem_name(pname,edgefile);


    if (debug) printf ("Debugging turned on\n");
    if (outfile) printf ("Output File: %s\n", outfile);
    fflush (stdout);


    rval = COLORread_dimacs (edgefile, &ncount, &ecount, &elist);
    COLORcheck_rval (rval, "COLORread_diamcs failed");

    tot_rtime = COLORcpu_time();

    if (cclasses_infile != (char*) NULL) {
       rval = COLORstable_read_stable_sets(&gcolors,&gcount,ncount,cclasses_infile,pname); 
       COLORcheck_rval(rval,"Failed in COLORstable_read_stable_sets");
    } else {
       rval = COLORgreedy (ncount, ecount, elist, &gcount, &gcolors);
       COLORcheck_rval (rval, "COLORgreedy failed");    

       /*     rval = COLORplot_graphviz(ncount,ecount,elist,0); */
       printf ("Greedy Colors: %d\n", gcount); fflush (stdout);
       for (i = 0; i < gcount; i++) {
          printf ("Color %d:", i);
          for (j = 0; j < gcolors[i].count; j++) {
             printf (" %d", gcolors[i].members[j]);
        }
          printf ("\n");
       }
       fflush (stdout);
    }
    
    gallocated = gcount;

    rval = COLORstable_initenv (&mwis_env,pname,write_mwis);
    COLORcheck_rval (rval, "COLORgreedy failed");

    rval = COLORlp_init (&lp, "colorme");
    COLORcheck_rval (rval, "COLORlp_init failed");

    for (i = 0; i < ncount; i++) {
       char sense = 'G';
        rval = COLORlp_addrow (lp, 0, (int *) NULL, (double *) NULL, sense,
                               1.0, (char*) NULL);
        COLORcheck_rval (rval, "COLORlp_addrow failed");
    }

    coef = (double *) malloc (ncount * sizeof (double));
    COLORcheck_NULL (coef, "out of memory for coef");
    for (i = 0; i < ncount; i++) coef[i] = 1.0;

    for (i = 0; i < gcount; i++) {
        rval = COLORlp_addcol (lp, gcolors[i].count, gcolors[i].members,
                               coef, 1.0, 0.0, 1.0, COLORlp_CONTINUOUS, (char*) NULL);
        if (rval) COLORlp_printerrorcode (rval);
        COLORcheck_rval (rval, "COLORlp_addcol failed");
    }

    rval = COLORlp_write (lp, "look.lp");
    COLORcheck_rval (rval, "COLORlp_write failed");

    pi = (double *) malloc (ncount * sizeof (double));
    COLORcheck_NULL (pi, "out of memory for pi");

    mwis_pi = (COLORNWT *) malloc (ncount * sizeof (COLORNWT));
    COLORcheck_NULL (mwis_pi, "out of memory for mwis_pi");

    
    do {
        ++ iterations;


        rval = COLORlp_optimize(lp);
        COLORcheck_rval (rval, "COLORlp_optimize failed");
       
        rval = COLORlp_pi (lp, pi);
        COLORcheck_rval (rval, "COLORlp_pi failed");

        make_pi_feasible(pi,gcolors,gcount);

        COLOR_double2COLORNWT(mwis_pi,&mwis_pi_scalef,
                              pi,ncount);

        print_objective(&lower_scaled_bound,mwis_pi,mwis_pi_scalef,ncount);

        {
            int set_i;

            rval = COLORstable_wrapper(&mwis_env,&newsets, &nnewsets, ncount, ecount,
                                       elist, mwis_pi,mwis_pi_scalef);
            COLORcheck_rval (rval, "COLORstable_gurobi failed");

            for (set_i = 0; set_i < nnewsets; ++set_i) {
                rval = COLORlp_addcol (lp, newsets[set_i].count,
                      newsets[set_i].members, coef, 1.0, 0.0, 1.0,
                      COLORlp_CONTINUOUS, NULL);
            }
            break_while_loop = (nnewsets == 0);

            concat_newsets(&gcolors,&gcount,&gallocated,
                           &newsets, &nnewsets);
        }
        
    } while ( (iterations < maxiterations) && !break_while_loop);

    if (iterations < maxiterations) {
       double incumbent;
       char   mps_fname[256];
       double dbl_lower_bound = COLORsafe_lower_dbl(lower_scaled_bound,mwis_pi_scalef);
       lower_bound = ceil(dbl_lower_bound);
       printf ("Found bound of %lld (%20.16g), greedy coloring %d (iterations = %d).\n", 
               (long long) lower_bound,dbl_lower_bound, gcount,iterations);

       tot_rtime = COLORcpu_time() - tot_rtime;
       printf("Computing initial lower bound took %f seconds.\n",tot_rtime);
       fflush(stdout);
       if (write_mwis) {
          sprintf(mps_fname,"%s.mwis.lp",pname);
          COLORstable_write_mps(mps_fname,ncount,ecount,elist,
                                mwis_pi,mwis_pi_scalef);

          sprintf(mps_fname,"%s.mwis.dimacs",pname);
          rval = COLORstable_write_dimacs(mps_fname,ncount,ecount,elist,
                                          mwis_pi,mwis_pi_scalef);

          sprintf(mps_fname,"%s.mwclq.dimacs",pname);
          rval = COLORstable_write_dimacs_clique(mps_fname,ncount,ecount,elist,
                                                 mwis_pi,mwis_pi_scalef);

          sprintf(mps_fname,"%s.graphviz.dot",pname);
          COLORplot_graphviz(mps_fname,ncount,ecount,elist,(int*) NULL); 
       }
       if (cclasses_outfile != (char*) NULL) {
          rval = COLORstable_write_stable_sets(gcolors,gcount,ncount,cclasses_outfile,pname); 
          COLORcheck_rval(rval,"Failed in COLORstable_write_stable_sets");
       }
       tot_rtime = COLORcpu_time();
       COLORlp_set_all_coltypes(lp,GRB_BINARY);
       COLORcheck_rval (rval, "COLORlp_set_all_coltypes");

       rval = COLORlp_optimize(lp);
       COLORcheck_rval (rval, "COLORlp_optimize failed");

       rval = COLORlp_objval (lp, &incumbent);
       COLORcheck_rval (rval, "COLORlp_objval failed");

       printf ("Found lower bound of %lld and upper bound of %g.\n", 
               (long long) lower_bound, incumbent);
       tot_rtime = COLORcpu_time() - tot_rtime;
       printf("Computing final upper bound took %f seconds.\n",tot_rtime);
       
       
    } else {
       printf ("Lower bound could not be found in %d iterations!\n", iterations);
    }
    
CLEANUP:

    if (lp) COLORlp_free (&lp);
    if (elist) free (elist);
    if (coef) free (coef);
    if (pi) free (pi);
    if (mwis_pi) free(mwis_pi);

    COLORfree_sets(&newsets,&nnewsets);
    COLORfree_sets(&gcolors,&gcount);
    COLORstable_freeenv(&mwis_env);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int rval = 0;

    while ((c = getopt (ac, av, "dmo:r:w:")) != EOF) {
        switch (c) {
        case 'd':
           /* each -d increases the verbosity by one.*/
           ++debug;
            break;
        case 'o':
            outfile = optarg;
            break;
        case 'm':
           write_mwis = 1;
           break;
        case 'r':
           cclasses_infile = optarg;
           break;
        case 'w':
           cclasses_outfile = optarg;
           break;
        default:
            usage (av[0]);
            rval = 1;  goto CLEANUP;
        }
    }

    if (ac <= optind) {
        rval = 1; goto CLEANUP;
    } else {
        edgefile = av[optind++];
    }

CLEANUP:

    if (rval) usage (av[0]);
    return rval;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage %s: [-see below-] edge_file\n", f);
    fprintf (stderr, "   -d    turn on debugging\n");
    fprintf (stderr, "   -o f  write coloring to file f\n");
    fprintf (stderr, "   -m    write final stable set and clique instances\n");
    fprintf (stderr, "   -r f  write initial stable sets from file f\n");
    fprintf (stderr, "   -w f  write stable sets from file f\n");
}

 static int get_problem_name(char* pname,const char* efname)
{
   int    rval = 0;
   int    len = 0;
   const char * fname = strrchr(efname,'/');
   char * lastdot = strrchr(efname,'.');
   if(!fname) {
      /* no slashes in efname.*/
      fname = efname;
   } else {
      fname++;
   }
   
   if (lastdot) {
      len = lastdot - fname + 1;
   } else {
      len = strlen(fname);
   }

   if (snprintf(pname,len,"%s",fname) < 0) {
      rval = 1;
   }
   printf("Extracted problem name %s\n",pname);

   return 0;
}

double COLORwall_time (void)
{
    return (double) time (0);
}

double COLORcpu_time (void)
{
    struct rusage ru;
    double t;

    getrusage (RUSAGE_SELF, &ru);

    t = ((double) ru.ru_utime.tv_sec) +
        ((double) ru.ru_utime.tv_usec) / 1000000.0;
    return t;
}

