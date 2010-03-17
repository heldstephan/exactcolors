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

#include "mwss.h"
#include "mwss_ext.h"


static int check_ncount(int ncount)
{
   int rval = 0;
   if (ncount > MAX_NODES) {
      fprintf(stderr,"SEWELL_init_graph ncount exceeds MAX_NODES.");
      rval = 1;
   }
   return rval;
}


static int SEWELL_init_graph(MWSSgraphpnt graph,
                             int ncount, int ecount,
                             const int elist[])
{
   int rval = 0;
   int      i, row, col;

   // Initialize the node names and degrees.
   graph->n_nodes = ncount;
   //MALLOC(node_list, n_nodes+1, tnode);
   for(i = 0; i <= graph->n_nodes; i++) {
      graph->node_list[i].name = i;
      graph->node_list[i].degree = 0;
   }

   // Initialize the adjacency matrix.
   for(row = 1; row <= graph->n_nodes; row++) {
      graph->adj[row][row] = 0;
      for(col = row + 1; col <= graph->n_nodes; col++) {
         graph->adj[row][col] = 0;
         graph->adj[col][row] = 0;
      }
   }

   for (i = 0; i < ecount; ++i) {
      row = elist[2*i]     + 1;
      col = elist[2*i + 1] + 1;
      graph->adj[row][col] = 1;
      graph->adj[col][row] = 1;
   }

   build_graph(graph);

   return rval;
}

extern
int SEWELL_optimize(int ** newset,
                    int   *nnewset,
                    int ncount, int ecount, const int elist[],
                    NWT nweights[]/* , */
/*                     NWT goal */)
{
   int rval = 0;
   wstable_info   info;
   int i;
   MWSSgraph      graph;
   MWSSdata       data;
   wstable_parameters parms;
   double density =  ((double) ecount) / ((double) (ncount * ( ncount - 1))) * 2.0;

   if(check_ncount(ncount)) {goto CLEANUP;}

   default_parameters(&parms);
   if (density < 0.2) {
      parms.clique_cover = 2;
   }

   SEWELL_init_graph(&graph,ncount,ecount,elist);
   for (i = 0; i < ncount; ++i) {
      graph.weight[i+1] = nweights[i];
   }


   initialize_max_wstable(&graph,&info);
   call_max_wstable(&graph,&data,&parms,&info);

   if (*newset) {free(*newset);}

   *nnewset = data.n_best;
   *newset  = (int*) malloc (data.n_best * sizeof(int));
   if (!*newset) {
      fprintf(stderr,"Failed to allocate *newset");
      rval = 1; goto CLEANUP;
   }

   for (i = 1; i <= data.n_best; ++i) {
      (*newset)[i - 1] = data.best_sol[i]->name - 1;
   }

 CLEANUP:
   free_graph(&graph);
   return rval;
}

