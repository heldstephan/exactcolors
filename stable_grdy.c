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

#include "color_defs.h"
#include "color.h"
#include "lp.h"
#include "graph.h"
#include "mwis.h"
#include "mwis_sewell/mwss_ext.h"

static char *graphfile = (char *) NULL;
static int debug = 0;
static double cpu_limit  = 0.0;
static int parseargs (int ac, char **av);
static int get_problem_name(char* pname,const char* efname);
static void usage (char *f);


static int parseargs (int ac, char **av)
{
    int c;
    int rval = 0;

    while ((c = getopt (ac, av, "dc:")) != EOF) {
       switch (c) {
       case 'd':
          debug = 1;
          break;
       case 'c':
          cpu_limit = atof(optarg);
          break;
       default:
          usage (av[0]);
          rval = 1;  goto CLEANUP;
       }
    }

    if (ac <= optind) {
        rval = 1; goto CLEANUP;
    } else {
        graphfile = av[optind++];
    }

CLEANUP:

    if (rval) usage (av[0]);
    return rval;
}
static int get_problem_name(char* pname,const char* efname)
{
    int rval = 0;
    int len = 0;
    const char *fname = strrchr(efname,'/');
    char *lastdot = strrchr(efname,'.');

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
        fprintf (stderr, "snprintf failed\n");
        rval = 1; goto CLEANUP;
    }
    printf("Extracted problem name %s\n",pname);  fflush (stdout);

CLEANUP:

   return rval;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage %s: [-see below-] graph_file\n", f);
    fprintf (stderr, "   -d       turn on debugging\n");
    fprintf (stderr, "   -c <double> number of seconds spend in B&B, negative values impose no limit\n");
}

int COLORdbg_lvl() {
   return debug;
}

int main (int ac, char **av)
{
    char pname[256] = "";
    int rval = 0;
    int ncount, ecount;
    int *elist = (int *) NULL;
    COLORNWT *nweights = (COLORNWT *) NULL;
    MWISls_env* ls_env   = (MWISls_env*) NULL;
    COLORset* newsets  = (COLORset*) NULL;
    int       nnewsets = 0;
    COLORNWT cutoff = 0;
    double szeit;
    int set_i, i;
    COLORNWT total_weight = 0;
    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    get_problem_name (pname, graphfile);

    rval = COLORread_dimacs (graphfile, &ncount, &ecount, &elist, &nweights);
    COLORcheck_rval (rval, "COLORread_dimacs failed");

    szeit = COLORutil_zeit ();


    rval = COLORstable_LS(&ls_env, &newsets, &nnewsets,
                          ncount, ecount, elist,
                          nweights, cutoff);
    COLORcheck_rval(rval,"COLORstable_LS failed.");

    for (set_i = 0; set_i < nnewsets; ++set_i) {
       total_weight = 0;
       for (i = 0; i < newsets[set_i].count; ++i) {
          int v = newsets[set_i].members[i];
          total_weight += nweights[v];
       }
       printf("Found stable set number %d with weight %d.\n",
              set_i + 1, total_weight);
    }

    printf("cpu_limit = %f\n", cpu_limit);
    if (cpu_limit > 0) {
       int*      bnb_newset  = (int*) NULL;
       int       bnb_nnewset = 0;
       COLORNWT  goal = COLORNWT_MAX;

       rval =  SEWELL_heur ( &bnb_newset, &bnb_nnewset,
                             ncount,ecount,
                             elist,  nweights,
                             total_weight,
                             goal,cpu_limit);

       if (rval == 0) {
          printf("B-&-B found maximum-weight stable set!\n");
       }

       if (rval == SEWELL_TIMEOUT) {
          printf("B-&-B reached cpu timeout!\n");
          rval = 0;
       }

       COLORcheck_rval(rval,"COLORstable_LS failed.");


       COLORNWT sewell_total_weight = 0;
          for (i = 0; i < bnb_nnewset; ++i) {
             int v = bnb_newset[i];
             sewell_total_weight += nweights[v];
          }
          if (sewell_total_weight > total_weight) {
             printf("Sewell:\n");
             printf("Found stable set number %d with weight %d.\n",
                    set_i + 1, sewell_total_weight);
          }
    }
    printf ("Running Time: %.2lf seconds\n", COLORutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:
    COLORstable_free_ls_env(&ls_env);

    COLOR_IFFREE(elist,int);
    COLOR_IFFREE(nweights,COLORNWT);
    COLOR_IFFREE(newsets,COLORset);

    return rval;
}
