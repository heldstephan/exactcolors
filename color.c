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
static int read_dimacs (char *f, int *pncount, int *pecount, int **pelist);


int COLORdbg_lvl() {
   return debug;
}

static void print_objective(double* pi,int ncount)
{
   int i;
   double obj = .0;
   for (i = 0; i < ncount;++i) {
      obj += pi[i];
   }
   printf("Current primal LP objective: %f.\n",obj);
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
         if  (pi[gcolors[c].members[i]] < 0) {
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
         printf("Decreased column sum of %5d from  %30.20f to  %30.20f\n",c,colsum,newcolsum);
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
    int rval = 0;
    int iterations = 0;
    int maxiterations = 1000000;
    int ncount = 0, ecount = 0, gcount, gallocated, i, j;
    int *elist = (int *) NULL;
    double *coef = (double *) NULL;
    double *pi = (double *) NULL;
    double lower_bound = 0;
    COLORset *gcolors = (COLORset *) NULL;
    COLORset *newsets = (COLORset*) NULL;   
    COLORlp *lp       = (COLORlp *) NULL;
    int       nnewsets = 0;
    int break_while_loop = 1;
    
    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    if (debug) printf ("Debugging turned on\n");
    if (outfile) printf ("Output File: %s\n", outfile);
    fflush (stdout);

    rval = read_dimacs (edgefile, &ncount, &ecount, &elist);
    COLORcheck_rval (rval, "read_diamcs failed");

    rval = COLORgreedy (ncount, ecount, elist, &gcount, &gcolors);
    COLORcheck_rval (rval, "COLORgreedy failed");

/*     rval = COLORplot_graphviz(ncount,ecount,elist,0); */
    gallocated = gcount;
    printf ("Greedy Colors: %d\n", gcount); fflush (stdout);
    for (i = 0; i < gcount; i++) {
        printf ("Color %d:", i);
        for (j = 0; j < gcolors[i].count; j++) {
            printf (" %d", gcolors[i].members[j]);
        }
        printf ("\n");
    }
    fflush (stdout);

    rval = COLORlp_init (&lp, "colorme");
    COLORcheck_rval (rval, "COLORlp_init failed");

    for (i = 0; i < ncount; i++) {
        rval = COLORlp_addrow (lp, 0, (int *) NULL, (double *) NULL, 'G',
                               1.0, NULL);
        COLORcheck_rval (rval, "COLORlp_addrow failed");
    }

    coef = (double *) malloc (ncount * sizeof (double));
    COLORcheck_NULL (coef, "out of memory for coef");
    for (i = 0; i < ncount; i++) coef[i] = 1.0;

    for (i = 0; i < gcount; i++) {
        rval = COLORlp_addcol (lp, gcolors[i].count, gcolors[i].members,
                          coef, 1.0, 0.0, 1.0, COLORlp_CONTINUOUS, NULL);
        if (rval) COLORlp_printerrorcode (rval);
        COLORcheck_rval (rval, "COLORlp_addcol failed");
    }

    rval = COLORlp_write (lp, "look.lp");
    COLORcheck_rval (rval, "COLORlp_write failed");

    pi = (double *) malloc (ncount * sizeof (double));
    COLORcheck_NULL (pi, "out of memory for pi");
    
    do {
        ++ iterations;


        rval = COLORlp_optimize(lp);
        COLORcheck_rval (rval, "COLORlp_optimize failed");
       
        rval = COLORlp_pi (lp, pi);
        COLORcheck_rval (rval, "COLORlp_pi failed");

        make_pi_feasible(pi,gcolors,gcount);

        print_objective(pi,ncount);

        {
            int set_i;

            rval = COLORstable_wrapper(&newsets, &nnewsets, ncount, ecount,
                                       elist, pi);
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
       rval = COLORlp_objval (lp, &lower_bound);
       COLORcheck_rval (rval, "COLORlp_objval failed");
    
       printf ("Found bound of %g (%g), greedy coloring %d (iterations = %d).\n", 
               ceil(lower_bound),lower_bound, gcount,iterations);

       COLORlp_set_all_coltypes(lp,GRB_BINARY);
       COLORcheck_rval (rval, "COLORlp_set_all_coltypes");

       rval = COLORlp_optimize(lp);
       COLORcheck_rval (rval, "COLORlp_optimize failed");

       rval = COLORlp_objval (lp, &incumbent);
       COLORcheck_rval (rval, "COLORlp_objval failed");

       printf ("Found lower bound of %g and upper bound of %g.\n", 
               ceil(lower_bound), incumbent);
       
       
    } else {
       printf ("Lower bound could not be found in %d iterations!\n", iterations);
    }

CLEANUP:

    if (lp) COLORlp_free (&lp);
    if (elist) free (elist);
    if (coef) free (coef);
    if (pi) free (pi);

    COLORfree_sets(&newsets,&nnewsets);
    COLORfree_sets(&gcolors,&gcount);
    
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
}

static int read_dimacs (char *f, int *pncount, int *pecount, int **pelist)
{
    int rval = 0;
    int ncount, ecount, icount = 0, haveprob = 0;
    int end0, end1;
    int *elist = (int *) NULL;
    char buf[256], *p;
    FILE *in = (FILE *) NULL;
    graph G;/* used to simplify graph.*/

    in = fopen (f, "r");
    if (!in) {
        fprintf (stderr, "Unable to open %s for input\n", f);
        rval = 1;  goto CLEANUP;
    }

    while (fgets (buf, 254, in) != (char *) NULL) {
        p = buf;
        if (p[0] == 'c') {
            printf ("Comment: %s", p+1);
        } else if (p[0] == 'p') {
           const char* delim = " \t\n";
           char* data = (char *) NULL;
            if (haveprob) {
                fprintf (stderr, "ERROR in Dimacs file -- two p lines\n");
                rval = 1;  goto CLEANUP;
            }
            haveprob = 1;
            data = strtok(p,delim); /* get 'p' */
            
            data = strtok(NULL,delim); /* get type */
            if ( strcmp(data,"edge") && strcmp(data,"edges") && strcmp(data,"col") ) {
               fprintf (stderr, "ERROR in Dimacs file -- not an edge file\n");
                rval = 1;  goto CLEANUP;
            }
            data = strtok(NULL,delim);
            sscanf (data, "%d", &ncount);
            data = strtok(NULL,delim);
            sscanf (data, "%d", &ecount);

            printf ("Number of Nodes: %d\n", ncount);
            printf ("Number of Edges: %d\n", ecount);
            elist = (int *) malloc (2 * ecount * sizeof (int));
            if (!elist) {
                fprintf (stderr, "out of memory for elist\n");
                rval = 1; goto CLEANUP;
            }
        } else if (p[0] == 'e') {
            if (!haveprob) {
                fprintf (stderr, "ERROR in Dimacs file -- e before p\n");
                rval = 1;  goto CLEANUP;
            }
            if (icount >= ecount) {
                fprintf (stderr, "ERROR in Dimacs file -- to many edges\n");
                rval = 1;  goto CLEANUP;
            }
            p++;
            sscanf (p, "%d %d", &end0, &end1);
/*             printf("Found edge %i %d %d\n", icount, end0,end1); */
            elist[2*icount] = end0-1;    /* Number nodes from 0, not 1 */
            elist[2*icount+1] = end1-1;
            icount++;
        }
    }

    rval = COLORadjgraph_build(&G, ncount,icount,elist);
    COLORcheck_rval(rval,"COLORadjgraph_build failed");                                     

    rval = COLORadjgraph_simplify(&G);
    COLORcheck_rval(rval,"COLORadjgraph_simplify failed");                                     

    COLORadjgraph_extract_edgelist(&icount, &elist,&G);
    COLORcheck_rval(rval,"COLORadjgraph_extract_edgelist");                                     

    *pncount = ncount;
    /* Some col-instances are buggy => reduce # edges to icount*/
    *pecount = icount; 
    *pelist = elist;

CLEANUP:
    COLORadjgraph_free(&G);
    if (rval) {
        if (elist) free (elist);
    }
    if (in) fclose (in);

    return rval;
}

