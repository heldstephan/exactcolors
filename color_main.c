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
#include <assert.h>

#include "color_private.h"

static char *edgefile = (char *) NULL;
static char *outfile = (char *) NULL;
static char *cclasses_infile = (char *) NULL;
static char *color_infile = (char *) NULL;

static int initial_upper_bound = INT_MAX;
static int parallel_branching  = 1;

static int get_problem_name(char* pname,const char* efname);

static void usage (char *f)
{
    fprintf (stderr, "Usage %s: [-see below-] edge_file\n", f);
    fprintf (stderr, "   -d    turn on debugging\n");
    fprintf (stderr, "   -o f  write coloring to file f\n");
    fprintf (stderr, "   -m    write final stable set and clique instances\n");
    fprintf (stderr, "   -r f  read initial stable sets from file f\n");
    fprintf (stderr, "   -w f  write stable sets from file f\n");
    fprintf (stderr, "   -c f  read initial coloring from file f\n");
    fprintf (stderr, "   -u int  initial upper bound f\n");
}


static int parseargs (int ac, char **av)
{
    int c;
    int rval = 0;
    int debug = COLORdbg_lvl();

    while ((c = getopt (ac, av, "dmo:r:w:c:u:")) != EOF) {
        switch (c) {
        case 'd':
           /* each -d increases the verbosity by one.*/
           ++debug;
           COLORset_dbg_lvl(debug);
           break;
        case 'o':
           outfile = optarg;
           break;
        case 'm':
           COLORset_write_mwis(1);
           break;
        case 'r':
           cclasses_infile = optarg;
           break;
        case 'w':
           COLORset_cclasses_outfile(optarg);
           break;
        case 'c':
           color_infile = optarg;
           break;
        case 'u':
           initial_upper_bound = atoi(optarg);
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


int main (int ac, char **av)
{
    int rval = 0;
    int i;
    double start_time;
    double tot_rtime;

    colordata  _colordata;
    colordata* cd = &_colordata;


    COLORset*  debugcolors = (COLORset*) NULL;
    int       ndebugcolors = 0;


    init_colordata(cd);

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    cd->upper_bound = initial_upper_bound;
    get_problem_name(cd->pname,edgefile);


    if (COLORdbg_lvl() > 1) printf ("Debugging turned on\n");
    if (outfile) printf ("Output File: %s\n", outfile);
    fflush (stdout);


    rval = COLORread_dimacs (edgefile, &(cd->ncount), &(cd->ecount),
                             &(cd->elist), (int **) NULL);
    COLORcheck_rval (rval, "COLORread_diamcs failed");
    cd->orig_node_ids = (int*) COLOR_SAFE_MALLOC(cd->ncount,int);
    COLORcheck_NULL(cd->orig_node_ids,"Failed to allocate cd->orig_node_ids");
    for (i = 0; i < cd->ncount; ++i) {cd->orig_node_ids[i] = i;}
    start_time = COLORcpu_time();

    if (cclasses_infile != (char*) NULL) {
       rval = COLORstable_read_stable_sets(&(cd->cclasses),&(cd->ccount),
                                           cd->ncount,cclasses_infile,cd->pname);
       COLORcheck_rval(rval,"Failed in COLORstable_read_stable_sets");
    } else {
       rval = COLORgreedy (cd->ncount, cd->ecount, cd->elist,
                           &(cd->ccount), &(cd->cclasses));
       COLORcheck_rval (rval, "COLORgreedycd failed");

       /*     rval = COLORplot_graphviz(ncount,ecount,elist,0); */
       printf ("Greedy Colors: %d\n", cd->ccount); fflush (stdout);
       print_colors(cd->cclasses,cd->ccount);
       COLORcopy_sets(&(cd->bestcolors),&(cd->nbestcolors),
                      cd->cclasses,cd->ccount);
       cd->upper_bound = cd->nbestcolors < cd->upper_bound ? cd->nbestcolors : cd->upper_bound;
    }

    if (color_infile != (char*) NULL) {
       rval = COLORstable_read_stable_sets(&debugcolors,&ndebugcolors,
                                           cd->ncount,color_infile,cd->pname);
       COLORcheck_rval(rval,"Failed in COLORstable_read_stable_sets");
       rval = COLORcheck_coloring(debugcolors,ndebugcolors,cd->ncount,
                           cd->ecount,cd->elist);
       COLORcheck_rval(rval,"Failed in COLORcheck_coloring");
       cd->debugcolors = debugcolors;
       cd->ndebugcolors = ndebugcolors;
       cd->opt_track = 1;
    }


    rval = compute_coloring(cd,parallel_branching);
    COLORcheck_rval(rval, "Failed to compute_coloring");

    if (cd->nbestcolors == cd->upper_bound) {
       printf ("Opt Colors: %d\n", cd->nbestcolors); fflush (stdout);
       print_colors(cd->bestcolors,cd->nbestcolors);
    } else {
       assert(cd->upper_bound == initial_upper_bound);
       printf("Lower bound reached predefined upper bound %d.\n",initial_upper_bound);
    }
    tot_rtime = COLORcpu_time() - start_time;
    printf("Computing coloring took %f seconds.\n",tot_rtime);


CLEANUP:
    if (debugcolors) free (debugcolors);
    free_colordata(cd);
    return rval;
}
