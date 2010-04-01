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
#include <float.h>
#include "mwss.h"

char     *prob_file;       // problem file

double   density = 0.0;    // -d option: density of graph (randomly generated graph or for printing)
int      gen = 0;          // -g option: randomly generate problem
int      lower_bound = 0;  // -l option: user supplied lower bound
int      n;                // -n option: # of nodes to randomly generate
double   seed = 3.1567;    // -s option: random seed (def = 3.1567)
int      test = 0;         // -t option: randomly generate and solve a set of test problems

/* Note: n, density, and seed can be supplied in order to randomly generate a problem,
         or they can be given in order to describe a problem which is read from a file
         (in which case these variables can be used for printing purposes.                  */


//_________________________________________________________________________________________________

int main(int ac, char **av)
{
   int            rval = 0;
   int            d;
   MWSSgraph      graph;
   MWSSdata       data;
   wstable_info   info;
   wstable_parameters parms;
   MWISNW         goal = MWISNW_MAX;     
   
   reset_pointers(&graph, &data, &info);

   default_parameters(&parms);

   parseargs (ac, av, &parms);

   if (test) {
      for(n = 10; n <= 100; n = n + 10) {
         for(d = 1; d <= 9; d = d + 2) {
            seed = 3.1567;
            rval = testprobs(&graph, &data, &parms, n, d / 10.0, &seed, &info);
            MWIScheck_rval(rval,"Failed in testprobs");
         }
      }
   } else {
      if (gen) {
         rval = rgrphgen(&graph, n, density, &seed);
         MWIScheck_rval(rval,"Failed in rgrphgen");
      } else {
         read_dimacs(&graph, prob_file);
      }
      if(parms.prn_info > 0) prn_graph(&graph);

      rval = initialize_max_wstable(&graph, &info);
      MWIScheck_rval(rval,"Failed in initialize_max_wstable");
      
/*       lower_bound = 6927366; */
      rval = call_max_wstable(&graph, &data, &parms, &info, goal, lower_bound);
      MWIScheck_rval(rval,"Failed in call_max_wstable");

   }

 CLEANUP:
   free_max_wstable(&graph,&data, &info);

   return rval;
}

//_________________________________________________________________________________________________

void parseargs(int ac, char **av, wstable_parameterspnt parms)
{
   int c, cnt;

   cnt = 0;
   while (++cnt < ac && av[cnt][0] == '-') {
      c = av[cnt][1];
      switch (c) {
         case 'c':
            parms->cpu_limit = atof(av[++cnt]);
            break;
         case 'd':
            density = atof(av[++cnt]);
            break;
         case 'g':
            gen = 1;
            break;
         case 'k':
            parms->clique_cover = atoi(av[++cnt]);
            break;
         case 'l':
            lower_bound = atoi(av[++cnt]);
            break;
         case 'n':
            n = atoi(av[++cnt]);
            break;
         case 'p':
            parms->prn_info = atoi(av[++cnt]);
            break;
         case 'r':
            parms->reorder = 1;
            break;
         case 's':
            seed = atof(av[++cnt]);
            break;
         case 't':
            test = 1;
            break;
         default:
            usage(*av);
            break;
      }
   }
   if (gen)
      prob_file = NULL;
   else {
      if (cnt >= ac) usage (*av);
      prob_file = av[cnt++];
   }
   if (cnt < ac) usage (*av);

   if (gen) {
      if (n <= 0) usage(*av);
      if (density < 0) usage(*av);
      if (seed == 0) usage(*av);
   }
}

//_________________________________________________________________________________________________

void usage (char *prog)
{
   fprintf (stderr, "Usage: %s probfile\n", prog);
   fprintf (stderr, "    -c: maximum cpu second (def=300)\n");
   fprintf (stderr, "    -d: density of graph (randomly generated graph or for printing)\n");
   fprintf (stderr, "    -g: randomly generate problem 1=yes 0=no (def = No)\n");
   fprintf (stderr, "    -k: clique cover 1 = maximal cliques  2 = maximum cliques (def = 1)\n");
   fprintf (stderr, "    -l: user supplied lower bound\n");
   fprintf (stderr, "    -n: # of nodes to randomly generate\n");
   fprintf (stderr, "    -p: controls level of printed information (def=0)\n");
   fprintf (stderr, "    -r: sort active nodes by ascending degree\n");
   fprintf (stderr, "    -s: seed for random number generation (def = 3.1567)\n");
   fprintf (stderr, "    -t: randomly generate and solve a set of test problems (def = No)\n");
   exit (1);
}
