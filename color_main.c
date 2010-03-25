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
    fprintf (stderr, "   -c f   read initial coloring from file f\n");
    fprintf (stderr, "   -p     start boss of parallel coloring\n");
    fprintf (stderr, "   -u int initial upper bound f\n");
    fprintf (stderr, "   -a     use B&B as coloring heuristic for upper bouns. f\n");
    fprintf (stderr, "   -l dbl cpu time limit for branching. f\n");

}


static int parseargs (int ac, char **av, COLORparms* parms)
{
    int c;
    int rval = 0;
    int debug = COLORdbg_lvl();

    while ((c = getopt (ac, av, "admpo:r:w:c:u:b:l:")) != EOF) {
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
        default:
           usage (av[0]);
           rval = 1;  goto CLEANUP;
        }
    }

    if (ac <= optind) {
        rval = 1; goto CLEANUP;
    } else {
       rval =  COLORparms_set_edgefile(parms,av[optind++]);
       COLORcheck_rval(rval,"Failed in COLORparms_set_edgefile");
    }

CLEANUP:

    if (rval) usage (av[0]);
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

   return 0;
}


static int print_hostinfo(void) {
    int   rval     = 0;

    char my_hostname[MAX_PNAME_LEN];
    pid_t my_pid   = getpid();

    rval = gethostname (my_hostname, MAX_PNAME_LEN - 1);
    COLORcheck_rval (rval, "gethostname failed");

    printf("Running color_worker on host %s with pid %lld.\n",
           my_hostname, (long long) my_pid);
 CLEANUP:
    return rval;
}


int main (int ac, char **av)
{
    int rval = 0;
    int i;
    double  start_time = COLORcpu_time();
    double tot_rtime;


    COLORproblem colorproblem;
    COLORparms* parms = &(colorproblem.parms);
    colordata*  cd = &(colorproblem.root_cd);


    COLORset*  debugcolors = (COLORset*) NULL;
    int        ndebugcolors = 0;

    print_hostinfo();

    COLORproblem_init(&colorproblem);
    cd->id = 0;
    colorproblem.ncolordata = 1;

    rval = parseargs (ac, av,parms);
    if (rval) goto CLEANUP;

    cd->upper_bound = parms->initial_upper_bound;
    get_problem_name(cd->pname,parms->edgefile);


    if (COLORdbg_lvl() > 1) printf ("Debugging turned on\n");
/*     if (parms->outfile)      printf ("Output File: %s\n", parms->outfile); */
    fflush (stdout);


    rval = COLORread_dimacs (parms->edgefile, &(cd->ncount), &(cd->ecount),
                             &(cd->elist), (int **) NULL);
    COLORcheck_rval (rval, "COLORread_diamcs failed");

    if (COLORget_backupdir()) {
       rval = recover_colordata(cd,&colorproblem);
    }
    if (cd->status == initialized) {
       cd->orig_node_ids = (int*) COLOR_SAFE_MALLOC(cd->ncount,int);
       COLORcheck_NULL(cd->orig_node_ids,"Failed to allocate cd->orig_node_ids");
       for (i = 0; i < cd->ncount; ++i) {cd->orig_node_ids[i] = i;}

       if (parms->cclasses_infile != (char*) NULL) {
          rval = COLORstable_read_stable_sets(&(cd->cclasses),&(cd->ccount),
                                              cd->ncount,parms->cclasses_infile,cd->pname);
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
          colorproblem.global_upper_bound = cd->upper_bound;
       }

       if (parms->color_infile != (char*) NULL) {
          rval = COLORstable_read_stable_sets(&debugcolors,&ndebugcolors,
                                              cd->ncount,parms->color_infile,cd->pname);
          COLORcheck_rval(rval,"Failed in COLORstable_read_stable_sets");
          rval = COLORcheck_coloring(debugcolors,ndebugcolors,cd->ncount,
                                     cd->ecount,cd->elist);
          COLORcheck_rval(rval,"Failed in COLORcheck_coloring");
          cd->debugcolors = debugcolors;
          cd->ndebugcolors = ndebugcolors;
          cd->opt_track = 1;
       }
    }

    rval = compute_coloring(&colorproblem);
    COLORcheck_rval(rval, "Failed to compute_coloring");

    if (cd->nbestcolors == cd->upper_bound) {
       printf ("Opt Colors: %d\n", cd->nbestcolors); fflush (stdout);
       print_colors(cd->bestcolors,cd->nbestcolors);
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

    if (debugcolors) free (debugcolors);

    return rval;
}
