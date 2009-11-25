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

typedef struct colordata {
   char  pname[256];

   /* The instance graph */
   int ncount;
   int ecount;
   int *elist;


   /* The column generation LP. */
   COLORlp * lp;
   double *coef;
   double *pi;
   COLORNWT  mwis_pi_scalef;
   COLORNWT *mwis_pi;
   COLORNWT  lower_bound;
   COLORNWT  lower_scaled_bound;

   /* The MWIS instances. */
   MWISenv*  mwis_env;
   COLORset *cclasses;
   int       ccount;
   int       gallocated;
   COLORset *newsets;
   int       nnewsets;
   
   
   int maxiterations;
   int retirementage;
} colordata;



int main (int ac, char **av);
static int parseargs (int ac, char **av);
static void usage (char *f);
static int get_problem_name(char* pname,const char* efname);
static int build_lp(colordata* cd);


static void init_colordata(colordata* cd)
{
   sprintf(cd->pname,"colorprob");
   cd->ncount = 0;
   cd->ecount = 0;
   cd->elist = (int *) NULL;

   cd->coef = (double *) NULL;
   cd->pi = (double *) NULL;
   cd->mwis_pi_scalef = 1.0;
   cd->mwis_pi = (COLORNWT *) NULL;
   cd->lower_bound = 0;
   cd->lower_scaled_bound = 0;
    
   cd->lp       = (COLORlp *) NULL;
   cd->mwis_env  = (MWISenv*) NULL;

   cd->newsets = (COLORset*) NULL;   
   cd->nnewsets = 0;

   cd->ccount = 0;
   cd->cclasses = (COLORset *) NULL;
   cd->gallocated = 0;      

   cd->maxiterations = 1000000;
   cd->retirementage = 1000000;
}
static void free_colordata(colordata* cd)
{
   if (cd->lp) COLORlp_free (&(cd->lp));
    if (cd->elist) free (cd->elist);
    if (cd->coef) free (cd->coef);
    if (cd->pi) free (cd->pi);
    if (cd->mwis_pi) free(cd->mwis_pi);

    COLORfree_sets(&(cd->newsets),&(cd->nnewsets));
    COLORfree_sets(&(cd->cclasses),&(cd->ccount));
    COLORstable_freeenv(&(cd->mwis_env));

}

int COLORdbg_lvl() {
   return debug;
}

static void print_objective(colordata* cd)
{
   int i;
   cd->lower_scaled_bound = .0;
   for (i = 0; i < cd->ncount;++i) {
      cd->lower_scaled_bound += (double) cd->mwis_pi[i];
   }

   printf("Current primal LP objective: %f (%lld / %lld).\n",
          COLORsafe_lower_dbl(cd->lower_scaled_bound,cd->mwis_pi_scalef),
                              (long long) cd->lower_scaled_bound,
                              (long long) cd->mwis_pi_scalef );
}

static void make_pi_feasible(colordata* cd)
{
   int c;
   int current_rounding = fegetround();
   fesetround(FE_UPWARD);

   for (c = 0; c < cd->ccount;++c) {
      int i;
      double colsum = .0;
      double newcolsum = .0;
      for (i = 0; i < cd->cclasses[c].count;++i) {
         if (signbit(cd->pi[cd->cclasses[c].members[i]])) {
            cd->pi[cd->cclasses[c].members[i]] = 0.0;
         }
         colsum += cd->pi[cd->cclasses[c].members[i]];
      }
      if (colsum > 1.0) {
         fesetround(FE_DOWNWARD);
         
         for (i = 0; i < cd->cclasses[c].count;++i) {
            cd->pi[cd->cclasses[c].members[i]] /= colsum;
            newcolsum += cd->pi[cd->cclasses[c].members[i]];
         }
         if (COLORdbg_lvl()> 1) {printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n",c,colsum,newcolsum);}
         fesetround(FE_UPWARD);
      }
      fesetround(current_rounding);
   }
}


static void reset_ages(COLORset* cclasses, int ccount) 
{
   int i;
   for (i = 0; i < ccount; ++i) {
      cclasses[i].age = 0;
   }
}

static int grow_ages(colordata* cd) 
{
   int rval = 0;
   int i;
   double* x = (double*) NULL;
   
   x = (double*) malloc(cd->ccount * sizeof(double));
   COLORcheck_NULL(x,"Failed to allocate x");

   rval = COLORlp_x(cd->lp,x);
   COLORcheck_rval(rval,"Failed in COLORlp_x");
   
   for (i = 0; i < cd->ccount; ++i) {
      if (x[i] < DBL_EPSILON) {
         ++(cd->cclasses[i].age);
      } else {
         cd->cclasses[i].age = 0;
      }
   }

 CLEANUP:
   if (x) free(x);
   return rval;
}

static int delete_old_colorclasses(colordata* cd) 
{
   int rval   = 0;
   int i;
   int numdel = 0;

   
   for (i = 0; i < cd->ccount; ++i) {
      if (cd->cclasses[i].age > cd->retirementage) {
         /* swap current class with last. */
         --(cd->ccount);
         memcpy(&(cd->cclasses[i]),&(cd->cclasses[cd->ccount]),sizeof(COLORset));
         /* Ensure that the formerly last class is considered:*/
         --i;
         numdel++;
      }
   }
   
   if (numdel) {
      printf("Deleted %d out of %d columns with age > %d. Rebuilding LP from scratch.\n",
             numdel, numdel + cd->ccount, cd->retirementage);
      COLORlp_free (&(cd->lp));

      rval = build_lp(cd);
      COLORcheck_rval(rval, "Failed in build_lp");
   }


 CLEANUP:
   
   return rval;
}

static int write_snapshot(colordata* cd, int add_timestamp) 
{
   int rval = 0;
   if (cclasses_outfile != (char*) NULL) {
      char   fname[256];

      if (add_timestamp) {
         char   timestr[256];
         time_t t;
         t = time(NULL);
         strftime(timestr, 256, "%Y%m%d%H%M", localtime(&t));
         sprintf(fname,"%s.%s",cclasses_outfile,timestr);
      } else {
         sprintf(fname,"%s",cclasses_outfile);
      }
      rval = COLORstable_write_stable_sets(cd->cclasses,cd->ccount,cd->ncount,
                                           fname,cd->pname); 
      COLORcheck_rval(rval,"Failed in COLORstable_write_stable_sets");
   }
 CLEANUP:
   return rval;
}


static void print_ages(colordata* cd) 
{
   int i;
   printf("AGES:");
   for (i = 0; i < cd->ccount; ++i) {
      printf(" %4d",cd->cclasses[i].age);
   }
   printf("\n");

}



static int concat_newsets(colordata* cd)
{
   int rval = 0;
   COLORset* tmpsets = (COLORset*) NULL;
   
   if (cd->nnewsets == 0) return rval;
   
   reset_ages(cd->newsets,cd->nnewsets) ;
      
   if (cd->ccount + cd->nnewsets > cd->gallocated) {
      /* Double size */
      tmpsets = malloc(2 * cd->gallocated * sizeof(COLORset));
      COLORcheck_NULL (tmpsets, "out of memory for tmpsets");
      memcpy(tmpsets,cd->cclasses, cd->ccount * sizeof(COLORset));
      free(cd->cclasses);
      cd->gallocated *= 2;
      cd->cclasses = tmpsets;
      tmpsets = NULL;
   }
   memcpy(cd->cclasses + cd->ccount, cd->newsets, cd->nnewsets * sizeof(COLORset));
   cd->ccount += cd->nnewsets;
 CLEANUP:
   if (rval) {
      if (cd->cclasses) free(cd->cclasses);
   }
   if (cd->newsets) free(cd->newsets);
   cd->newsets = (COLORset*) NULL;
   cd->nnewsets = 0;
   return rval;
}

static int build_lp(colordata* cd)
{
    int i;

    int rval = COLORlp_init (&(cd->lp), "colorme");
    COLORcheck_rval (rval, "COLORlp_init failed");

    for (i = 0; i < cd->ncount; i++) {
       char sense = 'G';
        rval = COLORlp_addrow (cd->lp, 0, (int *) NULL, (double *) NULL, sense,
                               1.0, (char*) NULL);
        COLORcheck_rval (rval, "COLORlp_addrow failed");
    }

    cd->coef = (double *) realloc (cd->coef,cd->ncount * sizeof (double));
    COLORcheck_NULL (cd->coef, "out of memory for coef");
    for (i = 0; i < cd->ncount; i++) cd->coef[i] = 1.0;

    for (i = 0; i < cd->ccount; i++) {
       rval = COLORlp_addcol (cd->lp, cd->cclasses[i].count, cd->cclasses[i].members,
                              cd->coef, 1.0, 0.0, 1.0, COLORlp_CONTINUOUS, (char*) NULL);
       if (rval) COLORlp_printerrorcode (rval);
       COLORcheck_rval (rval, "COLORlp_addcol failed");
    }

    rval = COLORlp_write (cd->lp, "look.lp");
    COLORcheck_rval (rval, "COLORlp_write failed");

    cd->pi = (double *) realloc (cd->pi,cd->ncount * sizeof (double));
    COLORcheck_NULL (cd->pi, "out of memory for pi");
 CLEANUP:
    if (rval) {
       COLORlp_free(&(cd->lp));
       if (cd->coef) free(cd->coef);
       if (cd->pi) free(cd->pi);
    }
    return rval;
}

int main (int ac, char **av)
{
    int rval = 0;
    int i,j;
    int iterations       = 0;
    int break_while_loop = 1;
    double start_time;
    double last_snapshot_time;

    double cur_time;

    double tot_rtime;
    
    colordata  _colordata;
    colordata* cd = &_colordata;

    init_colordata(cd);
    
    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    get_problem_name(cd->pname,edgefile);


    if (debug) printf ("Debugging turned on\n");
    if (outfile) printf ("Output File: %s\n", outfile);
    fflush (stdout);


    rval = COLORread_dimacs (edgefile,
                             &(cd->ncount), &(cd->ecount), &(cd->elist));
    COLORcheck_rval (rval, "COLORread_diamcs failed");

    start_time = COLORcpu_time();
    last_snapshot_time = start_time;

    if (cclasses_infile != (char*) NULL) {
       rval = COLORstable_read_stable_sets(&(cd->cclasses),&(cd->ccount),
                                           cd->ncount,cclasses_infile,cd->pname); 
       COLORcheck_rval(rval,"Failed in COLORstable_read_stable_sets");
    } else {
       rval = COLORgreedy (cd->ncount, cd->ecount, cd->elist,
                           &(cd->ccount), &(cd->cclasses));
       COLORcheck_rval (rval, "COLORgreedy failed");    
       
       /*     rval = COLORplot_graphviz(ncount,ecount,elist,0); */
       printf ("Greedy Colors: %d\n", cd->ccount); fflush (stdout);
       for (i = 0; i < cd->ccount; i++) {
          printf ("Color %d:", i);
          for (j = 0; j < cd->cclasses[i].count; j++) {
             printf (" %d", cd->cclasses[i].members[j]);
        }
          printf ("\n");
       }
       fflush (stdout);
    }
    
    reset_ages(cd->cclasses,cd->ccount) ;
    cd->gallocated = cd->ccount;

    rval = build_lp(cd);
    COLORcheck_rval (rval, "build_lp failed");

    rval = COLORstable_initenv (&(cd->mwis_env),cd->pname,write_mwis);
    COLORcheck_rval (rval, "COLORgreedy failed");

   
    cd->mwis_pi = (COLORNWT *) malloc (cd->ncount * sizeof (COLORNWT));
    COLORcheck_NULL (cd->mwis_pi, "out of memory for mwis_pi");
    
    cd->retirementage = cd->ncount + 1;
    do {
       ++iterations;

       if ( (iterations % (2 * cd->retirementage)) == 0) {
          delete_old_colorclasses(cd);
       }
       
       if (cclasses_outfile) {
          int add_timestamp = 1;
          cur_time = COLORcpu_time();
          if (cur_time - last_snapshot_time > 7200) {
             write_snapshot(cd,add_timestamp);
             last_snapshot_time = cur_time;
          }
       }

        rval = COLORlp_optimize(cd->lp);
        COLORcheck_rval (rval, "COLORlp_optimize failed");
       
        rval = grow_ages(cd);
        COLORcheck_rval (rval, "grow_ages failed");

        if (COLORdbg_lvl() > 1) {print_ages(cd);}

        rval = COLORlp_pi (cd->lp, cd->pi);
        COLORcheck_rval (rval, "COLORlp_pi failed");

        make_pi_feasible(cd);

        COLOR_double2COLORNWT(cd->mwis_pi,&(cd->mwis_pi_scalef),
                              cd->pi,cd->ncount);

        print_objective(cd);

        {
            int set_i;

            rval = COLORstable_wrapper(&(cd->mwis_env),&(cd->newsets), &(cd->nnewsets), 
                                       cd->ncount, cd->ecount,
                                       cd->elist, cd->mwis_pi,cd->mwis_pi_scalef);
            COLORcheck_rval (rval, "COLORstable_gurobi failed");
            
            for (set_i = 0; set_i < cd->nnewsets; ++set_i) {
                rval = COLORlp_addcol (cd->lp, cd->newsets[set_i].count,
                                       cd->newsets[set_i].members, 
                                       cd->coef, 1.0, 0.0, 1.0,
                                       COLORlp_CONTINUOUS, NULL);
            }
            break_while_loop = (cd->nnewsets == 0);

            concat_newsets(cd);
        }
        
    } while ( (iterations < cd->maxiterations) && !break_while_loop);

    if (iterations < cd->maxiterations) {
       double incumbent;
       char   mps_fname[256];
       double dbl_lower_bound = COLORsafe_lower_dbl(cd->lower_scaled_bound,cd->mwis_pi_scalef);
       cd->lower_bound = ceil(dbl_lower_bound);
       printf ("Found bound of %lld (%20.16g), greedy coloring %d (iterations = %d).\n", 
               (long long) cd->lower_bound,dbl_lower_bound, cd->ccount,iterations);

       tot_rtime = COLORcpu_time() - start_time;
       printf("Computing initial lower bound took %f seconds.\n",tot_rtime);
       fflush(stdout);
       if (write_mwis) {
          sprintf(mps_fname,"%s.mwis.lp",cd->pname);
          COLORstable_write_mps(mps_fname,cd->ncount,cd->ecount,cd->elist,
                                cd->mwis_pi,cd->mwis_pi_scalef);

          sprintf(mps_fname,"%s.mwis.dimacs",cd->pname);
          rval = COLORstable_write_dimacs(mps_fname,cd->ncount,cd->ecount,cd->elist,
                                          cd->mwis_pi,cd->mwis_pi_scalef);

          sprintf(mps_fname,"%s.mwclq.dimacs",cd->pname);
          rval = COLORstable_write_dimacs_clique(mps_fname,
                                                 cd->ncount,cd->ecount,cd->elist,
                                                 cd->mwis_pi,cd->mwis_pi_scalef);

          sprintf(mps_fname,"%s.graphviz.dot",cd->pname);
          COLORplot_graphviz(mps_fname,cd->ncount,cd->ecount,cd->elist,(int*) NULL); 
       }
       
       delete_old_colorclasses(cd);

       if (cclasses_outfile) {
          int add_timestamp = 0;
          write_snapshot(cd,add_timestamp);
       }

       tot_rtime = COLORcpu_time();
       COLORlp_set_all_coltypes(cd->lp,GRB_BINARY);
       COLORcheck_rval (rval, "COLORlp_set_all_coltypes");

       rval = COLORlp_optimize(cd->lp);
       COLORcheck_rval (rval, "COLORlp_optimize failed");

       rval = COLORlp_objval (cd->lp, &incumbent);
       COLORcheck_rval (rval, "COLORlp_objval failed");

       printf ("Found lower bound of %lld and upper bound of %g.\n", 
               (long long) cd->lower_bound, incumbent);
       tot_rtime = COLORcpu_time() - tot_rtime;
       printf("Computing final upper bound took %f seconds.\n",tot_rtime);
       
       
    } else {
       printf ("Lower bound could not be found in %d iterations!\n", iterations);
    }


    
CLEANUP:
    free_colordata(cd);
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

