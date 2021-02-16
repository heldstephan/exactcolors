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
#include <time.h>
#include <math.h>
#include <float.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>

#include "color_defs.h"
#include "color_parms.h"
#include "color_private.h"

static int get_problem_name(char* pname,const char* efname);

static void usage (char *f)
{
    fprintf (stderr, "Usage %s: [-see below-] edge_file\n", f);
    fprintf (stderr, "   -d     turn on debugging\n");
    fprintf (stderr, "   -b d   write intermediate solutions to directory d\n");
    fprintf (stderr, "   -o f   write coloring to file f\n");
    fprintf (stderr, "   -m     write final stable set and clique instances\n");
    fprintf (stderr, "   -r f   read initial stable sets from file f\n");
    fprintf (stderr, "   -w f   write stable sets to file f\n");
    fprintf (stderr, "   -c f   read initial coloring from file to track optimum path in the B&B-tree f\n");
    fprintf (stderr, "   -p     start boss of parallel coloring\n");
    fprintf (stderr, "   -u int initial upper bound\n");
    fprintf (stderr, "   -a     use B&B as coloring heuristic for upper bouns.\n");
    fprintf (stderr, "   -s int Branching strategy: 0 = none, 1 = minimum lower bound,"
             " 2 = DFS  (default), 3 = hybrid (2 followed by 1).\n");
    fprintf (stderr, "   -R int rounding style: 0 = neighbor (default), 1 = uniform, 2 = none\n");
    fprintf (stderr, "   -l dbl cpu time limit for branching.\n");

}


static int parseargs (int ac, char **av, COLORparms* parms)
{
    int c;
    int rval = 0;
    int debug = COLORdbg_lvl();

    while ((c = getopt (ac, av, "admpo:r:w:c:u:b:l:s:R:")) != EOF) {
        switch (c) {
        case 'd':
           /* each -d increases the verbosity by one.*/
           ++(debug);
           COLORset_dbg_lvl(debug);
           break;
        case 'o':
           rval = COLORparms_set_outfile(parms,optarg);
           COLORcheck_rval(rval,"Failed in COLORparms_set_outfile");
           break;
        case 'm':
           COLORparms_set_write_mwis(parms,1);
           break;
        case 'r':
           rval = COLORparms_set_cclasses_infile(parms,optarg);
           COLORcheck_rval(rval,"Failed in COLORparms_set_cclasses_infile");
           break;
        case 'w':
           rval = COLORparms_set_cclasses_outfile(parms,optarg);
           COLORcheck_rval(rval,"Failed in COLORparms_set_cclasses_outfile");
           break;
        case 'c':
           rval = COLORparms_set_color_infile(parms,optarg);
           COLORcheck_rval(rval,"Failed in COLORparms_set_color_infile");
           break;
        case 'u':
           rval = COLORparms_set_initial_upper_bound(parms,atoi(optarg));
           COLORcheck_rval(rval,"Failed in COLORparms_set_initial_upper_bound");
           break;
        case 'a':
           parms->upper_bounds_only  = 1;
           parms->branching_strategy = COLOR_hybrid_strategy;
           break;
        case 'p':
           rval = COLORparms_set_parallel(parms,1);
           COLORcheck_rval(rval,"Failed in COLORparms_set_initial_upper_bound");
           break;
        case 'b':
           rval = COLORparms_set_backupdir(parms,optarg);
           COLORcheck_rval(rval,"Failed in COLORparms_set_backupdir");
           break;
        case 'l':
           rval = COLORparms_set_branching_cpu_limit(parms,atof(optarg));
           COLORcheck_rval(rval,"Failed in COLORparms_set_branching_cpu_limit");
           break;
        case 's':
           rval = COLORparms_set_branching_strategy(parms,atoi(optarg));
           COLORcheck_rval(rval,"Failed in COLORparms_set_branching_strategy");
           break;
        case 'R':
           rval = COLORparms_set_rounding_strategy(parms,atoi(optarg));
           COLORcheck_rval(rval,"Failed in COLORparms_set_rounding_strategy");
           break;
        default:
          rval = 1;
           goto CLEANUP;
        }
    }

    if (ac <= optind) {
        rval = 1; goto CLEANUP;
    } else {
       rval =  COLORparms_set_edgefile(parms,av[optind++]);
       COLORcheck_rval(rval,"Failed in COLORparms_set_edgefile");
    }

CLEANUP:

    if (rval) {usage (av[0]);}
    return  (rval);
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

   return rval;
}

COLOR_MAYBE_UNUSED
static int quick_lower_bound(COLORset **newsets, int *nnewsets, int ncount,
                             int ecount, int *elist, int cutoff, int *pval) {
   int    rval = 0;
   int*   all_one_nweights = (int*) NULL;
   int    i;
   int    nrbranches = ncount/10 + 1;
   COLORset* cliques  = (COLORset*) NULL;
   int       ncliques;
   all_one_nweights = COLOR_SAFE_MALLOC(ncount, int);
   COLORcheck_NULL(all_one_nweights, "Failed to allocate all_one_nweights");

   for (i = 0; i < ncount; ++i) {all_one_nweights[i] = 1;}

   rval = COLORclique_ostergard(&cliques, &ncliques,
                                ncount, ecount, elist,
                                all_one_nweights,cutoff,pval,
                                nrbranches);
   COLORcheck_rval(rval,"Failed in COLORclique_ostergard.");


   /** From the found clique we initialize a coloring for the
       nodes in the cliques.
   */
   *nnewsets = cliques[ncliques-1].count;
   *newsets = COLOR_SAFE_MALLOC(*nnewsets,COLORset);
   COLORcheck_NULL(*newsets,"Failed to allocate *newsets");

   for (i = 0; i < cliques[ncliques-1].count; ++i) {
      (*newsets)[i].count = 1;
      (*newsets)[i].members = COLOR_SAFE_MALLOC(1,int);
      (*newsets)[i].members[0] = cliques[ncliques-1].members[i];
   }
 CLEANUP:
  COLORfree_sets(&cliques,&ncliques);
   COLOR_IFFREE(all_one_nweights,int);
   return rval;
}


int main (int ac, char **av)
{
    int rval = 0;
    double  start_time = COLORcpu_time();
    double tot_rtime;


    COLORproblem colorproblem;
    COLORparms* parms = &(colorproblem.parms);
    colordata*  cd = &(colorproblem.root_cd);

    int ncolors = 0;
    COLORset *colorclasses = (COLORset*) NULL;


    COLORset*  debugcolors = (COLORset*) NULL;

    rval = COLORprogram_header (ac,av);
    COLORcheck_rval(rval, "Failed in COLORprogram_header");

    COLORproblem_init(&colorproblem);
    cd->id = 0;
    colorproblem.ncolordata = 1;

    rval = parseargs (ac, av,parms);
    if (rval) goto CLEANUP;

    cd->upper_bound = parms->initial_upper_bound;
    get_problem_name(cd->pname,parms->edgefile);

    printf ("Rounding strategy: ");
     switch (parms->rounding_strategy) {
     case COLOR_neighbor_rounding:
        printf("neighbor\n");
        break;
     case COLOR_uniform_rounding:
        printf("uniform\n");
        break;
     case COLOR_no_rounding:
        printf("none\n");
        break;
     }

    if (COLORdbg_lvl() > 1) printf ("Debugging turned on\n");
    fflush (stdout);


    rval = COLORread_dimacs (parms->edgefile, &(cd->ncount), &(cd->ecount),
                             &(cd->elist), (int **) NULL);
    COLORcheck_rval (rval, "COLORread_diamcs failed");

    if (cd->upper_bound > cd->ncount) {cd->upper_bound = cd->ncount;}

    if (colorproblem.parms.backupdir) {
       recover_colordata(cd,&colorproblem);
    }

    rval = COLORexact_coloring(&colorproblem,&ncolors, &colorclasses);
    COLORcheck_rval(rval, "Failed to COLORexact_coloring");

    if (cd->nbestcolors == cd->upper_bound) {
       printf ("Opt Colors: %d\n", cd->nbestcolors); fflush (stdout);
       print_colors(colorclasses, ncolors);
    } else if (cd->lower_bound == parms->initial_upper_bound) {
       printf("Lower bound reached predefined upper bound %d.\n",
              parms->initial_upper_bound);
    } else {
       printf("Finished with LB %d and UB %d.\n",
              cd->lower_bound, cd->upper_bound);

    }
    tot_rtime = COLORcpu_time() - start_time;
    printf("Computing coloring took %f seconds.\n",tot_rtime);


CLEANUP:

    COLORproblem_free(&colorproblem);
    COLORfree_sets(&colorclasses, &ncolors);

    COLOR_IFFREE(debugcolors,COLORset);
    COLORlp_free_env();

    return rval;
}
