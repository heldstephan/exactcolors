#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <getopt.h>
#include <string.h>
#include <fenv.h>

#include "color.h"
#include "lp.h"
#include "graph.h"
#include "mwis.h"


static char *edgefile = (char *) NULL;
static char *outfile = (char *) NULL;
static int debug = 0;

int main (int ac, char **av);
static int parseargs (int ac, char **av);
static void usage (char *f);
static int get_problem_name(char* pname,const char* efname);

int main (int ac, char **av)
{
    char pname[256] = "";
    int rval = 0;

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    get_problem_name(pname,edgefile);


    if (debug) printf ("Debugging turned on\n");
    if (outfile) printf ("Output File: %s\n", outfile);
    fflush (stdout);

CLEANUP:

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int rval = 0;

    while ((c = getopt (ac, av, "do:")) != EOF) {
        switch (c) {
        case 'd':
            debug = 1;
            break;
        case 'o':
            outfile = optarg;
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

   return 0;
}

