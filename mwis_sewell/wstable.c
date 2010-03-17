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

/*
   1. This file was created on 12/11/09 by copying c:\sewell\research\stable\polar\concert\stable.cpp
   2. The purpose is to create an MWSS algorithm that Bill Cook can use to find improving stable sets
      in an IP graph coloring algorithm.
   3. Major modifications:
      a. Use as few global variables as possible.
      b. Use C intead of C++.
   4. Modified 3/5/10
      a. Eliminate static variables.
Previous History
   1. This file was created on 11/5/03 by modifying routines from stable.c, cliqcov.c, io.c, and
      sort.c from c:\sewell\research\stable\augment (which were copied and modified from
      c:\Linux\WASHU\research\independent\hybrid).
   2. Major modifications:
      a. Find an MWSS instead of an MSS.
      b. Use array positions 1,...,n_nodes instead of 0,...,n_nodes-1.
*/

#include <limits.h>

#include "mwss_ext.h"

#include "mwss.h"

//static   char     adj[MAX_NODES+1][MAX_NODES+1];  // adjacency matrix of the graph
                                                  // Note that it does not use first row or column
//static   int      n_nodes;                // number of nodes in the graph
//static   int      n_edges;                // number of edges in the graph
//static   tnode    node_list[MAX_NODES+1]; // vector of nodes of the graph
//static   nodepnt  *edge_list;             // space for list of edges
//static   nodepnt  **adj_last;             // points to last active node in adjacency list

//static   int      active_flag = 0;                    // used to test active status of a node
//static   nodepnt  cur_sol[MAX_NODES+1];               // Stores the current stable set during search
//static   nodepnt  act[MAX_NODES+1][MAX_NODES+1];      // Lists of active nodes at each level of search
//static   int      n_act[MAX_NODES+1];                 // n_active[d] = # of active nodes at depth d
//static   double   best_z;                             // Weight of best stable set found so far
//static   nodepnt  best_sol[MAX_NODES+1];              // Stores the best stable set found so far
//static   int      n_best;                             // # of nodes in best stable set found so far
//static   nodepnt  neighbors[MAX_NODES+1];             // work vector used by maximal_wclique
//static   int      count[MAX_NODES+1];                 // work vector used by ascending_distrib_sort
//static   nodepnt  list2[MAX_NODES+1];                 // work vector used by ascending_radix_sort

// Variables used by the Ostergard maximum clique/stable set routines.

//static   int   oster_cur_sol[MAX_NODES+1];            // Stores the current clique during search
//static   int   eligible[MAX_NODES+1][MAX_NODES+1];    // Lists of eligible nodes at each level of search
//static   int   found;                                 // found = 1 if a larger clique/stable set has been found
//static   int   n_elig[MAX_NODES+1];                   // n_eligible[d] = # of eligible nodes at depth d
//static   int   oster_best_z;                          // Number of nodes in largest clique/stable set found so far
//static   int   oster_best_sol[MAX_NODES+1];           // Stores the best clique/stable set found so far
//static   int   ub[MAX_NODES+1];                       // Used by oster_wclique.
//static   int   oster_neighbors[MAX_NODES+1];          // Work vector used by maximum_clique

//_________________________________________________________________________________________________

void initialize_max_wstable(MWSSgraphpnt graph, wstable_infopnt info)
/*
   1. This function should be called precisely once, prior to using max_wstable for the first time.
*/
{
   int   i;

   info->n_calls = 0;
   info->cpu = 0;
   info->clique_cover_cpu = 0;

   for(i = 1; i <= graph->n_nodes; i++) {
      //new_memory((void**) &node_list[i].adj_last, n_nodes+1, sizeof(nodepnt *));
      MALLOC(graph->node_list[i].adj_last, graph->n_nodes+1, nodepnt *);
      graph->node_list[i].adj_last[0] = graph->adj_last[i];
      graph->node_list[i].adj2 = graph->adj_last[i];
   }
}
//_________________________________________________________________________________________________
void default_parameters(wstable_parameterspnt parms)
{
   parms->clique_cover =   1; /* 1 = maximal cliques, 2 = maximum cliques*/
   parms->reorder      =   1;
   parms->cpu_limit    = 300;
   parms->prn_info     =   0;
}
//_________________________________________________________________________________________________

void call_max_wstable(MWSSgraphpnt graph, MWSSdatapnt data,
                      wstable_parameterspnt parameters,
                      wstable_infopnt info)
/*
   1. This routine sets up the data and calls max_wstable.
   2. Written 12/23/09.
   3. Modified 3/5/10.
      a. Accept a pointer to the graph.
      b. Eliminate weight as an input parameter.
      c. Accept a pointer to the data structures required for the search.
*/
{
   int                  i, n_best_stable, status;
   double               goal, lower_bound, z_best;
   nodepnt              *best_stable, *list;

   MALLOC(best_stable, graph->n_nodes + 1, nodepnt);
   MALLOC(list, graph->n_nodes + 1, nodepnt);
   for(i = 1; i <= graph->n_nodes; i++) list[i] = graph->node_list + i;
   goal = 0;
   for(i = 1; i <= graph->n_nodes; i++) goal += graph->weight[i];
   lower_bound = 0;
   //lower_bound = 8750000;
   status = max_wstable(graph, data, best_stable, &n_best_stable, &z_best, info,
                        parameters, list, graph->n_nodes, lower_bound, goal);

   free(best_stable);
   free(list);
}

//_________________________________________________________________________________________________

int max_wstable(MWSSgraphpnt graph, MWSSdatapnt data, nodepnt *best_stable, int *n_best_stable, double *z_best, wstable_infopnt info,
                wstable_parameterspnt parameters, nodepnt *list, int n_list, double lower_bound,
                double goal)
/*
   1. This routine searches for a weighted stable set in the graph induced
      by the nodes in list.  It only searches for stable sets with weight greater than goal.
      It initializes some data structures and then calls wstable to find the MWSS.
   2. See wstable for details about the search procedure.
   3. Input Variables
      a. best_stable stores the best stable set found so far.
      b. n_best_stable = # of nodes in best_stable.
      c. best_z = weight of best stable set found so far.
      d. info stores information about the search process.
      e. parameters contains the parameters for the search.
      f. list contains the nodes to be explored.
         Entries begin in position 1.
      g. n_list = # of nodes in list.
      h. lower_bound = lower bound for the problem.
      i. goal:  The search is terminated as soon as a stable set with weight >= goal is found.
         If you want to find an MWSS, then set goal equal to the sum of the node weights.
      j. weight[v] = weight of node v.
   4. Global Variables or static variables in this file
      a. adj is the adjacency matrix for the graph
      b. cur_sol stores the current stable set
      c. active is a list of the active nodes at each level of the search
      d. n_act[depth] = # of active nodes at depth d
      e. best_z = weight of best stable set found so far
                = lower_bound if a stable set with weight > lower_bound has not been found so far
      f. best_sol stores the best stable set found so far
      g. n_best = # of nodes in best stable set found so far
                = 0 if a stable set with weight > lower_bound has not been found so far
   5. Output Variables
      a. best_stable stores the best stable set found.
      b. n_best_stable = # of nodes in best_stable.
      c. z_best = weight of best stable set found.
      d. 1 is returned if a stable set with weight > lower_bound is found,
         0 o.w.
   6. Modified 3/5/10
      a. Accept a pointer to the graph.
      b. Eliminate weight as an input parameter.
      c. Accept a pointer to the data structures required for the search.
*/
{
   int      adj_last_offset, i, depth, n_active, n_out, status, v;
   double   cpu, greedy_weight, sum_weight;
   outnode  *out;
   osterdata   oster_data;

   depth = 1;
   n_active = n_list;
   sum_weight = 0;
   for(i = 1; i <= n_list; i++) {
      v = list[i]->name;
      assert((1 <= v) && (v <= graph->n_nodes));
      data->act[depth][i] = list[i];
      assert(graph->weight[v] >= 0);
      list[i]->weight = graph->weight[v];
      sum_weight += list[i]->weight;
      list[i]->adj_last[0] = list[i]->adj2;
   }
   data->n_act[depth] = n_active;
   for(i = 1; i <= graph->n_nodes; i++) graph->node_list[i].active = 0;
   graph->active_flag = 0;

   //reorder_nodes2(active, n_active);

   if(goal > sum_weight) goal = sum_weight;

   *n_best_stable = 0;
   *z_best = 0;
   data->best_z = lower_bound;
   data->n_best = 0;
   for(i = 0; i <= graph->n_nodes; i++) data->best_sol[i] = NULL;
   greedy_weight = greedy_wstable(graph, list, n_list, data->best_sol, &data->n_best);
   if(greedy_weight > data->best_z) {
      data->best_z = greedy_weight;
   }

   info->n_subproblems = 0;
   info->n_calls++;
   for(i = 0; i <= graph->n_nodes; i++) info->n_sub_depth[i] = 0;
   info->start_time = clock();

   adj_last_offset = 0;

   MALLOC(out, graph->n_nodes, outnode);
   n_out = 0;

   wstable(graph, data, &oster_data, adj_last_offset, 0, 1, goal, out, n_out, info, parameters);

   if(data->best_z > lower_bound) {
      *z_best = data->best_z;
      *n_best_stable = data->n_best;
      for(i = 1; i <= data->n_best; i++) best_stable[i] = data->best_sol[i];
      status = 1;
      //printf("best_z = %8.3f  greedy = %8.3f\n", best_z, greedy_weight);
      cpu = (double) (clock() - info->start_time) / CLOCKS_PER_SEC;
      if(parameters->prn_info) {
         printf("%3d %5d %8.0f %12d %10.3f %10.3f\n", 
                graph->n_nodes, graph->n_edges, data->best_z, 
                info->n_subproblems, info->clique_cover_cpu, cpu);
      }
   } else {
      status = 0;
   }

   info->cpu += (double) (clock() - info->start_time) / CLOCKS_PER_SEC;

   //printf("\nn = %d  m = %d  alpha = %d  n_sub = %d\n",
   //        n_nodes,n_edges, alpha, n_subproblems);
   //printf("-a = %d -k = %d -r = %d lb = %d ub = %d\n",
   //        adj_matrix,cover,reorder,lower_bound,upper_bound);

   if(parameters->prn_info) {
      printf("best_z = %6.3f\n", data->best_z);
      //prn_nodes(best_stable,n_best_stable);
      printf("n_subproblems = %d\n", info->n_subproblems);
      printf("depth n_sub_depth\n");
      for (i = 1; i <= data->n_best; i++)  printf("%5d %8d\n", i, info->n_sub_depth[i]);
   }

   free(out);

   return(status);
}

//_________________________________________________________________________________________________

void wstable(MWSSgraphpnt graph, MWSSdatapnt data, osterdatapnt oster_data, int adj_last_offset, double cur_z, int depth, double goal, outnode *out, int n_out,
           wstable_infopnt info, wstable_parameterspnt parameters)
//int      cover;          /* 1 to use clique_cover2,  2 to use clique_oc_cover */
//int      n_active;       /* number of node in active                           */
//int      *n_best_stable; /* number of node in the maximum stable set           */
//nodepnt  *active;        /* list of nodes to search for maximum stable set     *//
//nodepnt  *best_stable;   /* the nodes in a maximum stable set                  */
/*
   1. This subroutine searches for a stable set in the graph induced
      by the nodes in active[depth][*].  It only searches for stable sets with
      size greater than lower_bound.
   2. The search is terminated as soon as a stable set with size >= goal is found.
      If you want to find an MWSS, then set goal equal to the sum of the node weights.
   3. This routine assumes that the adjacency of two variables can be determined
      by adj[pntv->name][pntw->name], where pntv and pntw are pointers to the two nodes.
   4. This routine assumes the weights of the nodes are stored in pntv->weight.
      It also assumes that all the weights are nonnegative.
   5. Input variables
      a. adj_last_offset is used to determine when to stop reading the adjacency list.
         If the edges are reordered, then all the edges incident to active nodes
         in this subproblem will appear first in the adjacency list.
      b. cur_z = sum of weights of nodes in the current stable set.
      c. depth =  current depth in branch and bound search tree.
      d. goal:  The search is terminated as soon as a stable set with weight >= goal is found.
         If you want to find an MWSS, then set goal equal to the sum of the node weights.
      e. out is a list of nodes that have been fixed out of the current subproblem.
         Let v be a node in out.  Thus the stable set that we are searching for must contain
         neighbors of v whose node weights sum up to more than the weight of v.
         i.   This property can be used for pruning.
         ii.  This property can be used for choosing the branch nodes.
         iii. out does not necessarily contain every node that has been fixed out of the
              current solution.  Whenever enough of v's neighbors have been used, then
              v can be removed from out.
      f. n_out = # of nodes in out.
   6. Global Variables or static variables in this file
      a. adj is the adjacency matrix for the graph
      b. cur_sol stores the current clique
      c. active[depth][*] is a list of the active nodes for this subproblem
      d. n_act[depth] = # of active nodes for this subproblem
      e. best_z = weight of best stable set found so far
                = lower_bound if a stable set with weight > lower_bound has not been found so far
      f. best_sol stores the best stable set found so far
      g. n_best = # of nodes in best stable set found so far
                = 0 if a stable set with weight > lower_bound has not been found so far
   7. If a stable set with size greater than lower_bound is found, then
         a. The size of the largest stable set found is returned, and
         b. The largest stable set found is returned in best_stable.
      Otherwise,
         c. lower_bound is returned, and
         d. best_stable is not overwritten - it contains whatever it contained
            when this routine was called, which might be garbage.

   8. Warning:  The nodes in active[depth][*] may be sorted during this routine, so if the
      original order is important, it should be saved before calling.
   9. Warning: if out is used, then the subproblems must be
      generated by starting from the end of the active list.
  10. Modified 3/5/10.
      a. Accept a pointer to the graph.
      b. Accept a pointer to the data structures required for the search.
      c. Accept a pointer to the data structures required by oster_maximum_clique.
*/
{
   int      branch, j, n_active, n_branch_nodes, n_sub_active, n_sub_out, v;
   nodepnt  *active, branch_node;
   outnode  *sub_out;

   assert(cur_z >= 0.0);
   assert((1 <= depth) && (depth <= graph->n_nodes));
   //assert(sum_weight >= 0.0);
   n_active = data->n_act[depth];
   active = data->act[depth];
   assert((0 <= n_active) && (n_active <= graph->n_nodes));

   info->n_subproblems++;
   info->n_sub_depth[depth]++;

   // Check if a larger stable set has been found.

   if (cur_z > data->best_z) {
      data->best_z = cur_z;
      data->n_best = depth - 1;
      for(j = 1; j <= depth-1; j++) {
         v = data->cur_sol[j]->name;
         assert((1 <= v) && (v <= graph->n_nodes));
         data->best_sol[j] = data->cur_sol[j];
      }
      if(parameters->prn_info >= 1) {
         printf("\n Better stable set found.  depth = %3d best_z = %8.3f\n", depth, data->best_z);
         prn_nodes(data->best_sol, data->n_best);
      }
   }

   if (n_active == 0) {
      return;
   }

   reorder_nodes(graph, adj_last_offset, active, n_active);

   MALLOC(sub_out, n_active + n_out + 1, outnode);

   graph->active_flag++;
   for(j = 1; j <= n_active; j++) active[j]->active = graph->active_flag;
   if(reorder_out(graph, adj_last_offset/* , depth */, n_out, out) == 0) {
      //printf("p1 ");
      free(sub_out);
      return;
   }

   determine_branch_nodes(graph, data, oster_data, active, adj_last_offset, cur_z, n_active, &n_branch_nodes, out, n_out, info, parameters);

   //if((parameters->prn_info >= 2) || ((depth == 4) && (n_branch_nodes > 0))) {
   if(parameters->prn_info >= 2) {
      prn_stable(active, data->best_z, cur_z, active + n_active - n_branch_nodes, depth, n_active, n_branch_nodes, info, 1);
      //for(j = 1; j <= n_active; j++) printf("%3d ",active[j]->name); printf("\n");
   }


   for(j=n_active; (j>n_active-n_branch_nodes) && (data->best_z < goal); j--) {
   //for(j=n_active-n_branch_nodes+1; (j<=n_active) && (alpha<upper_bound); j++) {
      branch_node = active[j];
      stable_sub_problem(data, adj_last_offset, &branch, branch_node, cur_z, depth, j, &n_sub_active);
      if(branch) {
         new_out(adj_last_offset, branch_node, out, n_out, sub_out, &n_sub_out);

         wstable(graph, data, oster_data, adj_last_offset+1, cur_z + branch_node->weight, depth+1, goal, sub_out, n_sub_out, info, parameters);
      }
      out[++n_out].pntv = branch_node;
      out[n_out].w_in_cur_sol = 0;
      //if(prn_info >= depth) {
      //   printf("(%d %d %d %d %d %d) ",depth,branch_node->name+LOW_NODE,j,
      //           branch_node->key,n_sub_active,n_subproblems-n_save);
      //   fflush(stdout);
      //}
   }

   free(sub_out);
   return;
}

//_________________________________________________________________________________________________

void determine_branch_nodes(MWSSgraphpnt graph, MWSSdatapnt data, osterdatapnt oster_data, nodepnt *active, int adj_last_offset, double cur_z, int n_active,
                       int *n_branch_nodes, outnode *out, int n_out, wstable_infopnt info, wstable_parameterspnt parameters)
/*
   1. This routine determines the nodes on which to branch.
   2. Modified 3/5/10 to accept a pointer to the graph.
   3. Modified 3/12/10.
      a. Accept a pointer to the data structures required for the search.
      b. Accept a pointer to the data structures required by oster_maximum_clique.
*/
{
   int      i, n_active2;
   double   available_weight;
   nodepnt  *active2, *branch_nodes;

   MALLOC(branch_nodes, n_active+1, nodepnt);
   MALLOC(active2, n_active+1, nodepnt);

   *n_branch_nodes = 0;
   n_active2 = 0;
   available_weight = data->best_z - cur_z;

   weighted_clique_cover(graph, oster_data, active, active2, adj_last_offset, branch_nodes, n_active, &n_active2,
                         available_weight, n_branch_nodes, info, parameters);

   if(*n_branch_nodes > 1) branch_surplus(active, active2, branch_nodes, n_active, &n_active2, n_branch_nodes);

   //if(*n_branch_nodes >  1.0 * active[n_active]->key) {
   //   branch_neighbors(active, active2, adj_last_offset, branch_nodes, n_active, &n_active2, n_branch_nodes);
   //}

   if(n_active2 + *n_branch_nodes != n_active) {
      printf("Error: n_active2 + n_branch_nodes != n_active\n");
      exit(1);
   }

   branch_out(active, active2, branch_nodes, n_active, &n_active2, n_branch_nodes, out, n_out);

   ascending_distrib_sort(branch_nodes, graph->n_nodes, *n_branch_nodes);  // Don't sort branch nodes if using branch_neighbors
                                                                      // Must sort them inside of branch_neighbours

   for(i = 1; i <= *n_branch_nodes; i++) active2[++n_active2] = branch_nodes[i];
   for(i = 1; i <= n_active; i++) active[i] = active2[i];

   //cout << "n_branch_nodes = " << *n_branch_nodes << endl;
   //prn_nodes(active, n_active);

   free(branch_nodes);
   free(active2);
}

//_________________________________________________________________________________________________

void stable_sub_problem(MWSSdatapnt data, int adj_last_offset, int *branch, nodepnt branch_node, double cur_z,
                        int depth, int j, int *n_sub_active)
/*
   1. This routine creates a subproblem for the node to which branch_node points.
      It creates the list (sub_active) of nodes for the new subproblem.
   2. The branch flag is set to 1, unless the number of nodes in the subproblem
      is less than the lower bound for the subproblem.
   3. This routine uses an adjacency matrix to determine adjacencies.
   4. j must satisfy: branch_node = active[depth][j].
   5. Modified 3/12/10 to accept a pointer to the data structures required for the search.
*/
{
   char     *adjv;
   int      cnt, i;
   double   sub_weight;
   nodepnt  *active, *sub_active, pntw;

   assert(data->act[depth][j] == branch_node);
   //assert((1 <= branch_node->name) && (branch_node->name <= n_nodes));

   data->cur_sol[depth] = branch_node;

   active = data->act[depth];
   sub_active = data->act[depth + 1];
   adjv = branch_node->adjacent;

   *branch = 0;
   cnt = 0;
   *n_sub_active = 0;
   sub_weight = 0;
   for(i = 1; i < j; i++) {
      pntw = active[i];
      //assert((1 <= pntw->name) && (pntw->name <= n_nodes));
      if(!adjv[pntw->name]) {
         sub_active[++cnt] = pntw;
         pntw->adj_last[adj_last_offset+1] = pntw->adj_last[adj_last_offset];
         sub_weight += pntw->weight;
      }
   }
   *n_sub_active = cnt;
   data->n_act[depth+1] = cnt;

   if(cur_z + branch_node->weight + sub_weight > data->best_z) {
      *branch = 1;
   }
}

//_________________________________________________________________________________________________

void weighted_clique_cover(MWSSgraphpnt graph, osterdatapnt oster_data, nodepnt *active, nodepnt *active2, int adj_last_offset, nodepnt *branch_nodes,
              int n_active, int *n_active2, double available_weight, int *n_branch_nodes, wstable_infopnt info, wstable_parameterspnt parameters)
/*
   1. This routine finds a weighted clique cover for the nodes in active.
   2. The weight of the cliques in the cover is <= available_weight.
   3. The cover is created by processing the nodes in increasing order of their remaining weight.
      The remaining weight is the original weight minus the sum of the weights of the cliques
      in which it has already been included.
   4. Nodes do not explicitly receive a color, nor is the entire
      graph covered (a node whose remaining weight is > 0 is not covered).
      Thus, this routine does not develop an upper bound for the entire graph.
   5. Nodes that are not covered are placed in branch_nodes.
   6. Nodes that are covered are placed in active2.
   7. pntv->adj_last[adj_last_offset] is used to terminate reading the
      adjacency list of node v (node pointed to by pntv).  Thus, it must
      be set so that all of v's neighbors that are active are read.
   8. This routine initializes pntv->active to equal active_flag for all
      nodes in active.
   9. Created 11/7/03 by modifying clique_cover2 from c:\sewell\research\stable\augment.
  10. Modified 3/5/10 to accept a pointer to the graph.
  11. Modified 3/12/10 to accept a pointer to the data structures required by oster_maximum_clique.
*/
{
   int      i, n_act2;
   double   min_weight;
   nodepnt  min_node, pntv;
   clock_t  start_time;

   start_time = clock();

   graph->active_flag++;
   for(i = 1; i <= n_active; i++) {
      active[i]->active = graph->active_flag;
      active[i]->remaining_weight = active[i]->weight;
      active2[i] = active[i];
   }
   n_act2 = n_active;

   while(n_act2 > 0) {

      // Find the active node with the smallest remaining weight.
      // Remove covered nodes from active2.

      min_weight = 1.0E75;
      min_node = NULL;
      i = 1;
      while(i <= n_act2) {
         pntv = active2[i];

         if(pntv->remaining_weight < EPSILON) {
            pntv->active = 0;
            active2[i] = active2[n_act2];
            active2[n_act2] = NULL;
            n_act2--;
            continue;
         } else {
            i++;
         }

         if(pntv->remaining_weight < min_weight) {
            min_weight = pntv->remaining_weight;
            min_node = pntv;
         }
      }

      if(min_node == NULL) break;

      pntv = min_node;
      assert(pntv->remaining_weight > EPSILON);
      assert(pntv->active == graph->active_flag);
      pntv->active = 0;

      if(pntv->remaining_weight <= available_weight) {
         available_weight -= pntv->remaining_weight;
         if(parameters->clique_cover == 1) {
            maximal_wclique(graph, adj_last_offset, pntv);
         } else {
            maximum_clique(graph, oster_data, adj_last_offset, pntv);
         }
      } else {
         break;
      }
    }

   // Copy nodes into active2 and branch_nodes.

   n_act2 = 0;
   *n_branch_nodes = 0;
   for(i = 1; i <= n_active; i++) {
      pntv = active[i];
      if(pntv->remaining_weight < EPSILON) {
         active2[++n_act2] = pntv;
      } else {
         branch_nodes[++(*n_branch_nodes)] = pntv;
      }
   }
   *n_active2 = n_act2;

   info->clique_cover_cpu += (double) (clock() - start_time) / CLOCKS_PER_SEC;
}

//_________________________________________________________________________________________________

void maximal_wclique(MWSSgraphpnt graph, int adj_last_offset, nodepnt pntv)
/*
   1. This routine finds a maximal clique among the active neighbors of v.
   2. This routine assumes that pntw->active has already been set equal
      to active_flag for all desired nodes.
   3. Nodes are added to the clique in increasing order of pntw->remaining_weight, which must
      be initialized prior to calling this routine.
   4. This routine assumes that the adjancency matrix is available.
   5. This routine uses a global workvector, neighbors, which must
      have dimension at least equal to number of neighbors of v plus 1.
   6. Created 11/7/03 by modifying maximal_clique2 from c:\sewell\research\stable\augment.
   7. Modified 3/5/10 to accept a pointer to the graph.
   8. Modified 3/11/10.  Made neighbors a local variable.
*/
{
   double   min_weight, v_weight;
   nodepnt  min_node,pntw;
   char     *adjw;
   int      i,limit,n_neighbors;
   nodepnt  pntz,*ppnt,*ppnt_last;
   nodepnt  neighbors[MAX_NODES+1];

   min_node = NULL;
   n_neighbors = 0;
   v_weight = pntv->remaining_weight;
   //printf("(%10.0f: %3d", pntv->remaining_weight, pntv->name);

/* Initialize neighbors of v.  Find the first node to add to clique. */

   min_weight = 1.0E75;
   ppnt_last = pntv->adj_last[adj_last_offset];
   for (ppnt = pntv->adj; ppnt <= ppnt_last; ppnt++) {
      pntw = *ppnt;
      if (pntw->active == graph->active_flag) {
         neighbors[++n_neighbors] = pntw;
         if(pntw->remaining_weight < min_weight) {
            min_node = pntw;
            min_weight = pntw->remaining_weight;
         }
      }
   }

/* Repeatedly remove w's neighbors, and find the next node
   to add to the clique. */

   pntw = min_node;
   while (pntw != NULL ) {
      pntw->remaining_weight -= v_weight;
      //printf(" %3d", pntw->name);
      adjw = pntw->adjacent;
      min_node = NULL;
      min_weight = 1.0E75;
      limit = n_neighbors;
      n_neighbors = 0;
      for (i = 1; i <= limit; i++) {
         pntz = neighbors[i];
         if (adjw[pntz->name]) {
            neighbors[++n_neighbors] = pntz;
            if(pntz->remaining_weight < min_weight) {
               min_node = pntz;
               min_weight = pntz->remaining_weight;
            }
         }
      }
      pntw = min_node;
   }
   pntv->remaining_weight = 0;
   //printf(")\n");
}

//_________________________________________________________________________________________________

void maximum_clique(MWSSgraphpnt graph, osterdatapnt oster_data, int adj_last_offset, nodepnt pntv)
/*
   1. This routine finds a maximum clique among the active neighbors of v.
   2. This routine assumes that pntw->active has already been set equal
      to active_flag for all desired nodes.
   3. oster_maximum_clique is used to find a maximum clique.
   4. This routine assumes that the adjancency matrix is available.
   5. This routine uses a global workvector, oster_neighbors, which must
      have dimension at least equal to number of neighbors of v plus 1.
   6. Created 1/15/10 by modifying maximal_wclique.
   7. Modified 3/5/10 to accept a pointer to the graph.
   8. Modified 3/12/10 to accept a pointer to the data structures required by oster_maximum_clique.
*/
{
   double   v_weight;
   nodepnt  pntw;
   int      i, n_neighbors;
   nodepnt  *ppnt, *ppnt_last;

   n_neighbors = 0;
   v_weight = pntv->remaining_weight;
   //printf("(%10.0f: %3d", pntv->remaining_weight, pntv->name);

   // Initialize neighbors of v.

   ppnt_last = pntv->adj_last[adj_last_offset];
   for (ppnt = pntv->adj; ppnt <= ppnt_last; ppnt++) {
      pntw = *ppnt;
      if (pntw->active == graph->active_flag) {
         oster_data->oster_neighbors[++n_neighbors] = pntw->name;
      }
   }

   // Call oster_maximum_clique to find a maximum clique among the neighbors of v.

   oster_data->oster_best_z = 0;
   if(n_neighbors > 0) oster_maximum_clique(graph, oster_data, oster_data->oster_neighbors, n_neighbors, 1);

   // Reduce the weight of all the nodes in the maximum clique.

   for(i = 1; i <= oster_data->oster_best_z; i++) {
      graph->node_list[oster_data->oster_best_sol[i]].remaining_weight -= v_weight;
   }
   pntv->remaining_weight = 0;
   //printf(")\n");
}

//_________________________________________________________________________________________________

void oster_clique(MWSSgraphpnt graph, osterdatapnt oster_data, int cur_z, int depth, char clique_or_stable)
/*
   1. This is a recursive procedure to find a maximum clique/stable set.
   2. It uses Ostergard's maximum clique algorithm.
   3. Input Variables
      a. cur_z = number of nodes in the current solution
      b. depth = depth in branch and bound search tree
      c. clique_or_stable = '1' for maximum weight clique, '0' for mwss
   4. Global Variables or static variables in this file
      a. adj is the adjacency matrix for the graph
      b. oster_cur_sol stores the current clique
      c. eligible[depth,*] is a list of the eligible nodes at each level of the search
      d. n_elig[depth] = # of eligible nodes at depth d
      e. found = 1 if a larger clique/stable set has been found
      f. oster_best_z = size of best clique/stable set found so far
      g. oster_best_sol stores the best clique/stable set found so far
      h. ub[v_k] = upper bound on size of maxium clique/stable set in graph induced by nodes
                   v_1, ... , v_k, where these are the nodes that have been
                   examined thus far in the B&B tree.
   5. Created 1/13/10 by modifying oster_wclique in c:\sewell\research\stable\polar\concert\cliqstab.cpp
   6. Modified 3/5/10 to accept a pointer to the graph.
   7. Modified 3/12/10 to accept a pointer to the data structures required by oster_maximum_clique.
*/
{
   char  *adjv;
   int   cnt, depth1, i, j, v, w;

   assert(cur_z >= 0);
   assert((1 <= depth) && (depth <= graph->n_nodes));
   assert((clique_or_stable == 0) || (clique_or_stable == 1));
   assert((0 <= oster_data->n_elig[depth]) && (oster_data->n_elig[depth] <= graph->n_nodes));

   if (oster_data->n_elig[depth] == 0) {
      if (cur_z > oster_data->oster_best_z) {
         oster_data->oster_best_z = cur_z;
         for(i = 1; i <= depth-1; i++) oster_data->oster_best_sol[i] = oster_data->oster_cur_sol[i];
         oster_data->found = 1;
      }
      return;
   }

   depth1 = depth + 1;
   for(i = oster_data->n_elig[depth]; i > 0; i--) {
      if (cur_z + i <= oster_data->oster_best_z) return;
      v = oster_data->eligible[depth][i];
      assert((1 <= v) && (v <= graph->n_nodes));
      adjv = graph->node_list[v].adjacent;
      oster_data->oster_cur_sol[depth] = v;
      assert((0 <= oster_data->ub[v]) && (oster_data->ub[v] <= oster_data->oster_best_z));
      if (cur_z + oster_data->ub[v] <= oster_data->oster_best_z) return;
      cnt = 0;
      for(j = 1; j <= i-1; j++) {
         w = oster_data->eligible[depth][j];
         assert((1 <= w) && (w <= graph->n_nodes));
         if (adjv[w] == clique_or_stable) {
            oster_data->eligible[depth1][++cnt] = w;
         }
      }
      oster_data->n_elig[depth1] = cnt;
      oster_clique(graph, oster_data, cur_z + 1, depth1, clique_or_stable);
   }
}

//_________________________________________________________________________________________________

int oster_maximum_clique(MWSSgraphpnt graph, osterdatapnt oster_data, int list[], int n_list, char clique_or_stable)
/*
   1. This procedure uses oster_clique to find a maximum clique or stable set.
   2. Input Variables
      a. list contains the nodes that are to be explored.
         Entries begin in position 1.
      b. n_list = # of variables in list.
      c. clique_or_stable = '1' for maximum weight clique, '0' for mwss.
   3. Global Variables - see oster_clique
   4. Output Variables
      a. The size of the maximum clique or stable set is returned.
      b. The maximum clique or stable set is not returned, but it is available in oster_best_sol.
   5. Created 1/15/10 by modifying maximum_wclique in c:\sewell\research\stable\polar\concert\cliqstab.cpp
   6. Modified 3/5/10 to accept a pointer to the graph.
   7. Modified 3/12/10 to accept a pointer to the data structures required by oster_maximum_clique.
*/
{
   char  *adjv;
   int   cnt, i, j, v, w;

   // Set the ub for each node in list.

   for (i = 1; i <= n_list; i++) {
      v = list[i];
      oster_data->eligible[1][i] = v;
      //ub[v] = 0;
   }

   oster_data->oster_best_z = 0;
   for(i = 1; i <= n_list; i++) {
      oster_data->found = 0;
      v = oster_data->eligible[1][i];
      adjv = graph->node_list[v].adjacent;
      //n_elig[1] = i;
      oster_data->oster_cur_sol[1] = v;
      cnt = 0;
      for(j = 1; j <= i-1; j++) {
         w = oster_data->eligible[1][j];
         //if (adj[v][w] == clique_or_stable) {
         if (adjv[w] == clique_or_stable) {
            oster_data->eligible[2][++cnt] = w;
         }
      }
      oster_data->n_elig[2] = cnt;
      oster_clique(graph, oster_data, 1, 2, clique_or_stable);
      oster_data->ub[v] = oster_data->oster_best_z;
   }
   assert(check_cliq(graph, oster_data, clique_or_stable, oster_data->oster_best_z));

   return(oster_data->oster_best_z);
}

//_________________________________________________________________________________________________

int check_cliq(MWSSgraphpnt graph, osterdatapnt oster_data, char clique_or_stable, int n_best)
/*
   1. This routine checks the clique (or stable set) found by oster_wlique.
      It does not verify that	it is a maximum clique (or stable set).
	2. 1 is returned if everything is OK, o.w. 0 is returned.
   3. Modified 3/5/10 to accept a pointer to the graph.
   4. Modified 3/12/10 to accept a pointer to the data structures required by oster_maximum_clique.
*/
{
   int    i, j, v;

   assert((clique_or_stable == 0) || (clique_or_stable == 1));
   assert((1 <= n_best) && (n_best <= graph->n_nodes));

	for(i = 1; i <= n_best; i++) {
      v = oster_data->oster_best_sol[i];
      assert((1 <= v) && (v <= graph->n_nodes));
	   for(j = i+1; j <= n_best; j++) {
         if(graph->adj[v][oster_data->oster_best_sol[j]] != clique_or_stable) {
				printf("oster_clique did not find a clique\n");
				return(0);
			}
		}
	}
   return(1);
}

//_________________________________________________________________________________________________

void branch_surplus(nodepnt *active, nodepnt *active2, nodepnt *branch_nodes, int n_active,
                    int *n_active2, int *n_branch_nodes)
/*
   1. This routine checks if any of the active nodes have surplus <= 0.  If so, one such node is
      made the branch node.  If no such node is found, then active2 and branch_nodes are not modified.
   2. This routine assumes the current surplus of each node is stored ->surplus.
*/
{
   int      i, index;

   index = -1;
   for(i = 1; i <= n_active; i++) {
      if(active[i]->surplus <= 0) {
         index = i;
         break;
      }
   }

   if(index > 0) {
      *n_active2 = 0;
      for(i = 1; i <= n_active; i++) {
         if(i != index) {
            //assert((1 <= active[i]->name) && (active[i]->name <= n_nodes));
            active2[++(*n_active2)] = active[i];
         }
      }
      *n_branch_nodes = 1;
      branch_nodes[*n_branch_nodes] = active[index];
   }


}

//_________________________________________________________________________________________________

void branch_neighbors(nodepnt *active, nodepnt *active2, /* int adj_last_offset, */ nodepnt *branch_nodes,
                      int n_active, int *n_active2, int *n_branch_nodes)
/*
   1. This routine chooses a node v and uses it and its neighbors as the branch nodes.
   2. Warning: this routine assumes the current degree of each node is stored ->key.
*/
{
   char     *adjv;
   int      i, index;
   nodepnt  pntv, pntw;

   index = n_active;

   pntv = active[index];
   adjv = pntv->adjacent;
   *n_active2 = 0;
   *n_branch_nodes = 0;
   for(i = 1; i <= n_active; i++) {
      if(i != index) {
         pntw = active[i];
         //assert((1 <= pntw->name) && (pntw->name <= n_nodes));
         if(adjv[pntw->name]) {
            branch_nodes[++(*n_branch_nodes)] = pntw;
         } else {
            active2[++(*n_active2)] = pntw;
         }
      }
   }

   //ascending_distrib_sort(branch_nodes,n_nodes,*n_branch_nodes);
   branch_nodes[++(*n_branch_nodes)] = pntv;
}

//_________________________________________________________________________________________________

void reorder_nodes(MWSSgraphpnt graph, int adj_last_offset, nodepnt *list, int n_list)
/*
   1. This routine sorts the nodes in list in ascending order of degree in
      the graph induced by list.
   2. The degree of each node is calculated.  During the calculation, the
      adjacency list of each node in list is shuffled so that all of its
      neighbors in list appear first.
   3. list[i]->adj_last[adj_last_offset] is used to terminate reading
      the adjacency list of node i (node pointed to by list[i]).
      Thus, list[i]->adj_last[adj_last_offset] must be set so that all of
      i's neighbors in list will be read.
   4. After the degree of node i has been computed,
      list[i]->adj_last[adj_last_offset] is set to point to the last neighbor
      of i that is in list.
   5. The surplus of each node is computed and stored in ->surplus.
   6. Created 11/7/03 by modifying reorder_nodes from c:\sewell\research\stable\augment.
   7. Modified 3/12/10 to accept a pointer to the graph.
*/
{
   int      active_flag, i;
/*    int      ascending_key(); */
   double   surplus;
   nodepnt  pnt,pntu,*ppnt1,*ppnt2,*ppnt_last;

   graph->active_flag++;
   active_flag = graph->active_flag;
   //for(i = 1; i <= n_nodes; i++) node_list[i].active = 0;
   for(i = 1; i <= n_list; i++) list[i]->active = active_flag;

   for(i = n_list; i > 0; i--) {
      pntu = list[i];
      surplus = -list[i]->weight;
      ppnt_last = pntu->adj_last[adj_last_offset];
      for(ppnt1 = (ppnt2 = pntu->adj)-1; ppnt2 <= ppnt_last; ppnt2++) {
         if((*ppnt2)->active == active_flag) {
            surplus += (*ppnt2)->weight;
            pnt = *(++ppnt1);
            *ppnt1 = *ppnt2;
            *ppnt2 = pnt;
         }
      }
      pntu->adj_last[adj_last_offset] = ppnt1;
      pntu->key = (ppnt1 - pntu->adj) + 1;
      pntu->surplus = surplus;
   }

   //prn_subgraph(adj_last_offset, list, n_list);
   //assert(check_subgraph(adj_last_offset, list, n_list));

   //ascending_radix_sort(list,n_list);
}

void reorder_edges(MWSSgraphpnt graph, int adj_last_offset, nodepnt *list, int n_list)
/*
   1. The adjacency list of each node in list is shuffled so that all of its
      neighbors in list appear first.
   2. The full adjacency list is used.
   3.
   4. After the degree of node i has been computed,
      list[i]->adj_last[adj_last_offset] is set to point to the last neighbor
      of i that is in list.
   5. The surplus of each node is computed and stored in ->surplus.
   6. Created 11/19/03 by modifying reorder_nodes from c:\sewell\research\stable\polar\concert\stable.cpp.
      a. Read the full adjacency list rather than using adj_last[adj_last_offset] to terminate
         reading the adjacency list.
      b. Do not store the degree in key.
      c. Do not compute the surplus.
      d. Do not sort the nodes.
      e. Set list[i]->tabu_adj to point to the last neighbor of i that is in the list.
   7. Modified 3/12/10 to accept a pointer to the graph.
*/
{
   int      active_flag, i;
   nodepnt  pnt,pntu,*ppnt1,*ppnt2;

   graph->active_flag++;
   active_flag = graph->active_flag;
   for(i = 1; i <= n_list; i++) list[i]->active = active_flag;

   for(i = n_list; i > 0; i--) {
      pntu = list[i];
      for(ppnt1 = (ppnt2 = pntu->adj)-1; *ppnt2 != NULL; ppnt2++) {
         if((*ppnt2)->active == active_flag) {
            pnt = *(++ppnt1);
            *ppnt1 = *ppnt2;
            *ppnt2 = pnt;
         }
      }
      pntu->adj_last[adj_last_offset] = ppnt1;
      //pntu->tabu_adj = ppnt1;
   }

   //prn_subgraph(adj_last_offset, list, n_list);
   //assert(check_subgraph(adj_last_offset, list, n_list));
}

//_________________________________________________________________________________________________

void branch_out(nodepnt *active, nodepnt *active2, nodepnt *branch_nodes, int n_active,
                int *n_active2, int *n_branch_nodes, outnode *out, int n_out)
/*
   1. This routine searches out to find the branch nodes.
   2. Warning: this routine assumes the current degree of each node in out is is stored ->key.
   3. Created 11/13/03 by modifying branch_one_two from c:\sewell\research\stable\augment\cliqcov.c
*/
{
   char     *adjv;
   int      i, min_deg;
   nodepnt  min_node, pntv, pntw;

   min_deg = INT_MAX;
   min_node = NULL;

   for(i = 1; i <= n_out; i++) {
      pntv = out[i].pntv;
      if(pntv->key < min_deg) {
         min_deg = pntv->key;
         min_node = pntv;
      }
   }

   //if((min_node != NULL) && (*n_branch_nodes > min_deg) && (min_deg < 2)) {
   if((min_node != NULL) && (*n_branch_nodes > min_deg)) {
      adjv = min_node->adjacent;
      *n_active2 = 0;
      *n_branch_nodes = 0;
      for(i = 1; i <= n_active; i++) {
         pntw = active[i];
         if(adjv[pntw->name]) {
            branch_nodes[++(*n_branch_nodes)] = pntw;
         } else {
            active2[++(*n_active2)] = pntw;
         }
      }
      assert(min_node->key == *n_branch_nodes);

      //ascending_distrib_sort(branch_nodes, n_nodes, *n_branch_nodes);
   }
}

//_________________________________________________________________________________________________

int reorder_out(MWSSgraphpnt graph, int adj_last_offset/* , int depth */, int n_out, outnode *out)
/*
   1. The adjacency list of each node in out is shuffled so that all of its
      active neighbors appear first.
   2. This routine assumes that pntv->active has been set to active_flag for all
      active nodes and is less than active_flag for all inactive nodes.
   3. out[i]->adj_last[adj_last_offset] is used to terminate reading
      the adjacency list of node i (node pointed to by out[i]).
      Thus, out[i]->adj_last[adj_last_offset] must be set so that all of
      i's active neighbors will be read.
   4. The degree of each node is computed and one[i]->adj_last[adj_last_offset]
      is set to point to the last active neighbor of node i.
   5. For each node v in out, this routine computes the sum of the weights of v's active neighbors
      plus the sum of the weights of v's neighbors in cur_sol.  If this sum is less than or equal
      to the weight of v, then the current subproblem is infeasible and 0 is returned immediately.
      Otherwise 1 is returned.
   6. Created 11/12/03 by modifying reorder_one from c:\sewell\research\stable\augment\sort.c
   7. Modified 3/12/10 to accept a pointer to the graph.
*/
{
   int      active_flag, branch, i;
   double   sum_weight;
   nodepnt  pntu, pntw, *ppnt1, *ppnt2, *ppnt_last;

   active_flag = graph->active_flag;
   branch = 1;
   for(i = 1; i <= n_out; i++) {
      pntu = out[i].pntv;
      sum_weight = 0;
      ppnt_last = pntu->adj_last[adj_last_offset];
      for(ppnt1 = (ppnt2 = pntu->adj)-1; ppnt2 <= ppnt_last; ppnt2++) {
         pntw = *ppnt2;
         if(pntw->active == active_flag) {
            *ppnt2 = *(++ppnt1);
            *ppnt1 = pntw;
            sum_weight += pntw->weight;
         }
      }
      pntu->adj_last[adj_last_offset] = ppnt1;
      pntu->key = (ppnt1 - pntu->adj) + 1;

      if(pntu->weight >= out[i].w_in_cur_sol + sum_weight) return(0);
   }
   return(branch);
}

//_________________________________________________________________________________________________

void new_out(int adj_last_offset, nodepnt branch_node, outnode *out, int n_out, outnode *sub_out, int *n_sub_out)
/*
   1. This routine places the nodes of out that satisfy the following condition in sub_out.
      a. The sum of the weights of v's neighbors in cur_sol is <= v's weight.
   2. Created 11/12/03 by modifying new_one from c:\sewell\research\stable\augment\stable.c
*/
{
   char     *adjv;
   int      i;
   double   branch_weight;
   nodepnt  pnt;

   adjv = branch_node->adjacent;
   branch_weight = branch_node->weight;

   *n_sub_out = 0;
   for(i = 1; i <= n_out; i++) {
      pnt = out[i].pntv;
      if(!adjv[pnt->name]) {
         assert(out[i].w_in_cur_sol <= pnt->weight);
         sub_out[++(*n_sub_out)] = out[i];
         pnt->adj_last[adj_last_offset+1] = pnt->adj_last[adj_last_offset];
      } else {
         if(out[i].w_in_cur_sol + branch_weight <= pnt->weight) {
            sub_out[++(*n_sub_out)] = out[i];
            sub_out[*n_sub_out].w_in_cur_sol += branch_weight;
            pnt->adj_last[adj_last_offset+1] = pnt->adj_last[adj_last_offset];
         }
      }
   }
}

//_________________________________________________________________________________________________

void ascending_distrib_sort(nodepnt *list, int m, int n)
/*
   1. This routine sorts list in ascending order of key.
   2. It assumes that 0 <= key < m.
   3. n = # of items in list.
   4. This routine uses two global workvectors, count and list2.
      The dimension of count must be at least m, and the dimension
      of list2 must be at least n.
   5. Created 11/7/03 by modifying ascending_distrib_sort from c:\sewell\research\stable\augment.
   6. Modified 3/11/10.  Made count and list2 local variables instead of global variables.
*/
{
   int      i,i2,sum;
   int      count[MAX_NODES+1];
   nodepnt  list2[MAX_NODES+1];


   for(i = 0; i < m; i++) count[i] = 0;
   for(i = 1; i <= n; i++) {
      list2[i] = list[i];
      count[list[i]->key]++;
   }

   sum = 0;
   i = 0;
   while(sum < n) {
      sum += count[i];
      count[i++] = sum;
   }

   for(i = n; i > 0; i--) {
      i2 = count[list2[i]->key]--;
      list[i2] = list2[i];
   }
}

//_________________________________________________________________________________________________

double greedy_wstable(MWSSgraphpnt graph, nodepnt *list, int n_list, nodepnt *stable_set, int *n_stable_set)
/*
   1. This routine finds a stable set in the graph induced by the nodes in list.
      It repeatedly finds a node of minimum surplus, places that node in stable_set,
      and removes that node's neighbors from list.
   2. The weight of the stable set is returned.
   3. This routine replaces an earlier version of greedy_stable.  This one
      uses a different algorithm to update the degrees after each iteration.
      The new running time is O(m), where m is the number of edges in the graph
      induced by list.
   4. Uses list[i]->surplus to store surplus.
   5. Created 10/29/99 by modifying greedy_stable from
      linux\washu\research\independent\greedy.c
      a. Conform to ansi C.
      b. Read through entire adjacency list instead of using
         list[i]->adj_last[adj_last_offset].
      c. list2 and neighbors are static arrays instead of allocating each call.
      d. Don't use reorder_nodes (which set the active status,calculated the
         degrees, shuffled the adjacency lists, and sorted the nodes.
   6. Created 10/19/03 by modifying greedy_stable from c:\sewell\research\stable\gstabu\greedy.c
      a. list2 and neighbors are allocated in this routine.
   7. Created 11/10/03 by modifying greedy_stable from c:\sewell\research\stable\augment\augment.c.
   8. Modified 3/12/10 to accept a pointer to the graph.
*/
{
   int      i, limit, n_active, n_list2, n_neighbors;
   double   alpha_w, min_surplus, v_weight;
   nodepnt  *list2, min_node, *neighbors, *ppnt, pntv, pntw;

   MALLOC(list2, MAX_NODES+1, nodepnt);
   MALLOC(neighbors, MAX_NODES+1, nodepnt);

   n_active = n_list;
   n_list2 = n_list;
   *n_stable_set = 0;
   alpha_w = 0;
   if(n_list > 0) {
      graph->active_flag++;
      //for(i = 1; i <= n_nodes; i++) {
      //   node_list[i].active = 0;
      //   node_list[i].surplus = -node_list[i].weight;
      //}
      for(i = 1; i <= n_list; i++) {
         list[i]->active = graph->active_flag;
         list[i]->surplus = -list[i]->weight;
         list2[i] = list[i];
      }

      // Calculate surplus of each node (store in ->surplus)

      for(i = n_list; i > 0; i--) {
         pntv = list[i];
	      for (ppnt = pntv->adj; (pntw = *ppnt) != NULL; ppnt++) {
            if( pntw->active == graph->active_flag ) {
               pntv->surplus += pntw->weight;
            }
         }
      }


      while(n_active > 0) {

         // Find the node with minimum surplus and reduce list2 to active nodes.

         min_surplus = 1.0E75;
         min_node = NULL;
         limit = n_list2;
         n_list2 = 0;
         for(i = 1; i <= limit; i++) {
            pntv = list2[i];
            if(pntv->active == graph->active_flag) {
               if(pntv->surplus < min_surplus) {
                  min_node = pntv;
                  min_surplus = pntv->surplus;
               }
               list2[++n_list2] = pntv;
            }
         }
         n_active = n_list2;

         if(min_node != NULL) {
            stable_set[++(*n_stable_set)] = min_node;
            alpha_w += min_node->weight;
            min_node->active = 0;
            n_active--;

/*          inactivate min_node's neighbors and place them in neighbors */

            n_neighbors = 0;
            for(ppnt = min_node->adj; (pntv = *ppnt) != NULL; ppnt++) {
               if( pntv->active == graph->active_flag) {
                  pntv->active = 0;
                  n_active--;
                  neighbors[++n_neighbors] = pntv;
               }
            }

/*          for each neighbor, v, of min_node, decrement the degree of all
            of v's neighbors.  */

            for(i = 1; i <= n_neighbors; i++) {
               pntv = neighbors[i];
               v_weight = pntv->weight;
               for(ppnt = pntv->adj; (pntw = *ppnt) != NULL; ppnt++) {
                  if(pntw->active == graph->active_flag) {
                     pntw->surplus -= v_weight;
                  }
               }
            }
         }
      }

   }

   free(list2);
   free(neighbors);

   return(alpha_w);
}

//_________________________________________________________________________________________________

void build_graph(MWSSgraphpnt graph)
{
   int      i, j;
   nodepnt  *p;

   graph->n_edges = 0;
   for(i = 1; i < graph->n_nodes; i++) {
      for(j = i+1; j <= graph->n_nodes; j++) {
         if(graph->adj[i][j] == 1) {
            graph->n_edges++;
         }
      }
   }

   //MALLOC(node_list, n_nodes+1, tnode);
   MALLOC(graph->edge_list, (2*graph->n_edges)+graph->n_nodes, nodepnt);
   MALLOC(graph->adj_last, graph->n_nodes+1, nodepnt*);
   graph->node_list[0].adjacent = NULL;
   for(i = 1; i <= graph->n_nodes; i++) {
      graph->node_list[i].adjacent = graph->adj[i];
   }

   for(i = 0; i <= graph->n_nodes; i++) {
      graph->node_list[i].name = i;
      graph->node_list[i].degree = 0;
      graph->node_list[i].adjv = 0;
      graph->node_list[i].adj2 = NULL;
   }
   graph->node_list[0].name = -1;

   for(i = 0, p = graph->edge_list; i < ((2*graph->n_edges)+graph->n_nodes); i++,p++) {
      *p = NULL;
   }

   for(i = 1; i < graph->n_nodes; i++) {
      for(j = i+1; j <= graph->n_nodes; j++) {
         if(graph->adj[i][j] == 1) {
            assert(i != j);
            graph->node_list[i].degree++;
            graph->node_list[j].degree++;
         }
      }
   }

   p = graph->edge_list;
   for(i = 1; i <= graph->n_nodes; i++) {
      graph->node_list[i].adj = p;
      p += (graph->node_list[i].degree + 1);
   }

   for(i = 1; i < graph->n_nodes; i++) {
      for(j = i+1; j <= graph->n_nodes; j++) {
         if(graph->adj[i][j] == 1) {
            *(graph->node_list[i].adj) = &(graph->node_list[j]);
            if(graph->node_list[i].adj >= p) {
               fprintf(stderr, "out of bounds\n");
            }
            graph->node_list[i].adj++;
            *(graph->node_list[j].adj) = &(graph->node_list[i]);
            if(graph->node_list[j].adj >= p) {
               fprintf(stderr, "out of bounds\n");
            }
            graph->node_list[j].adj++;
         }
      }
   }

   p = graph->edge_list;
   for(i = 1; i <= graph->n_nodes; i++) {
      graph->node_list[i].adj = p;
      graph->adj_last[i] = p + graph->node_list[i].degree - 1;
      p += (graph->node_list[i].degree + 1);
   }
}

//_________________________________________________________________________________________________

int check_graph(MWSSgraphpnt graph)
{
   int      cnt, i, j;
   nodepnt  *ppnt, pntv;


   for(i = 1; i <= graph->n_nodes; i++) {
      for(j = 1; j <= graph->n_nodes; j++) {
         if(graph->node_list[i].adjacent[j] != graph->adj[i][j]) return(0);
         if(graph->node_list[j].adjacent[i] != graph->adj[i][j]) return(0);
      }
   }

   for(i = 1; i <= graph->n_nodes; i++) {
      cnt = 0;
      for(ppnt = graph->node_list[i].adj; *ppnt != NULL; ppnt++) {
         pntv = *ppnt;
         if(graph->adj[i][pntv->name] != 1) return(0);
         cnt++;
      }
      if(cnt != graph->node_list[i].degree) return(0);
   }

   return(1);
}

//_________________________________________________________________________________________________

int check_subgraph(MWSSgraphpnt graph, int adj_last_offset, nodepnt *list, int n_list)
{
   int      i, j, v, w;
   nodepnt  *ppnt, *ppnt_last, pntv, pntw;

   // Check that every node on each adjacency list is active.

   graph->active_flag++;
   for(i = 1; i <= n_list; i++) {
      list[i]->active = graph->active_flag;
   }

   for(i = 1; i <= n_list; i++) {
      pntv = list[i];
      ppnt_last = pntv->adj_last[adj_last_offset];
      for (ppnt = pntv->adj; ppnt <= ppnt_last; ppnt++) {
         pntw = *ppnt;
         if(pntw->active != graph->active_flag) {
            return(0);
         }
      }
   }

   // Check that the appropriate edges are on the adjacency lists for every pair of
   // nodes in the list that are adjacent.

   for(i = 1; i <= n_list; i++) {

      // Mark all the nodes in list that are on the adjacency list of v.

      graph->active_flag++;
      pntv = list[i];
      v = pntv->name;
      ppnt_last = pntv->adj_last[adj_last_offset];
      for (ppnt = pntv->adj; ppnt <= ppnt_last; ppnt++) {
         pntw = *ppnt;
         pntw->active = graph->active_flag;
      }

      // Check that every node on list that is adjacent to v has been marked.

      for(j = 1; j <= n_list; j++) {
         pntw = list[j];
         w = pntw->name;
         if(graph->adj[v][w] == 1) {
            if(pntw->active != graph->active_flag) {
               return(0);
            }
         }
      }
   }

   return(1);
}

//_________________________________________________________________________________________________

void prn_graph(MWSSgraphpnt graph)
{
   int      cnt,i;
   nodepnt  *pnt;
   nodepnt  pnt2;

   printf("\n");
   printf("%5d %5d\n", graph->n_nodes, graph->n_edges);
   for(i = 1; i <= graph->n_nodes; i++) {
      printf("%5d%5d",i, graph->node_list[i].degree);
      cnt = 0;
      pnt = graph->node_list[i].adj;
      while(*pnt != NULL) {
         cnt++;
         pnt2 = *pnt;
         printf("%5d%s",pnt2->name, (*(++pnt) != NULL) && (cnt % 14) == 0 ? "\n          ":"");
      }
      //if((cnt % 14) != 0) printf("\n");
      printf("\n");
   }
   for(i = 1; i <= graph->n_nodes; i++) {
      printf("%3d %10.0f\n", i, graph->weight[i]);
   }
}

//_________________________________________________________________________________________________

void prn_subgraph(int adj_last_offset, nodepnt *list, int n_list)
{
   int      cnt,i;
   nodepnt  *ppnt, *ppnt_last, pntv, pntw;

   for(i = 1; i <= n_list; i++) {
      pntv = list[i];
      printf("%4d:", pntv->name);
      cnt = 0;
      ppnt_last = pntv->adj_last[adj_last_offset];
      for (ppnt = pntv->adj; ppnt <= ppnt_last; ppnt++) {
         cnt++;
         pntw = *ppnt;
         printf(" %4d", pntw->name);
         if((cnt % 14) == 0) {
            printf("\n     ");
         }
      }
      printf("\n");
   }
   printf("\n");
}

//_________________________________________________________________________________________________

void prn_nodes(nodepnt *list, int n)
{
   int      i;

   for (i = 1; i <= n; i++) {
      printf("%5d%s",list[i]->name, (i % 16) == 0 ? "\n":"");
   }
   if ( (i % 16) != 0 )  printf("\n");
}

//_________________________________________________________________________________________________

void prn_node_weights(nodepnt *list, int n)
{
   int      i;

   for (i = 1; i <= n; i++) {
      printf(" %5.3f", list[i]->weight);
      if((i % 16) == 0) printf("\n");
   }
   if ( (i % 16) != 0 ) printf("\n");
}

//_________________________________________________________________________________________________

void prn_stable(nodepnt *active, double best_z, double cur_z, nodepnt *branch_nodes, int depth,
                int n_active, int n_branch_nodes, wstable_infopnt info, int prn_level)
{
   int      i;
   double   cpu, ratio, ratio2, sum_w, sum_uncovered;

   sum_w = 0;
   for(i = 1; i <= n_active; i++) sum_w += active[i]->weight;
   sum_uncovered = 0;
   for(i = 1; i <= n_branch_nodes; i++) sum_uncovered += branch_nodes[i]->weight;
   ratio = (sum_w - sum_uncovered) / (best_z - cur_z);
   ratio2 = sum_w / (best_z - cur_z);
   cpu = (double) (clock() - info->start_time) / CLOCKS_PER_SEC;
   if(prn_level > 1) {
      printf("depth = %3d n_sub = %10d  cpu = %8.2f cur_z = %10.0f best_z = %10.0f  sum_w = %10.0f n_active = %3d sum_covered = %10.0f n_branch_nodes = %3d  ratio = %5.2f  ratio2 = %5.2f\n",
              depth, info->n_subproblems, cpu, cur_z, best_z, sum_w, n_active, sum_w - sum_uncovered, n_branch_nodes, ratio, ratio2);
   } else {
      printf("%3d %10d  %8.2f %10.0f %10.0f  %10.0f %3d %10.0f %3d %5.2f %5.2f \n", depth, info->n_subproblems, cpu, cur_z, best_z, sum_w, n_active, sum_w - sum_uncovered, n_branch_nodes, ratio, ratio2);
   }
   if(prn_level >= 3) {
      printf("active\n");
      prn_nodes(active, n_active - n_branch_nodes);
      prn_node_weights(active, n_active - n_branch_nodes);
      printf("branch_nodes \n");
      prn_nodes(branch_nodes, n_branch_nodes);
      prn_node_weights(branch_nodes, n_branch_nodes);
   }
}


/*************************************************************************************************/

/*
   1. The following functions randomly generate problems.
   2. They were created 3/18/08 by modifying the corresponding routines from
      c:\sewell\research\stable\dbfs\dbfs\dbfs.cpp.
*/

void rgrphgen(MWSSgraphpnt graph, int n, double density, double* dseed)
/*
   1. This function generate a random graph, where the probability of
      an edge being included equals density.
   2. Input Variables
      a. n = number of nodes in the graph.
      b. density = The probability that an edge is created between two nodes.
      c. dseed = seed used for random number generation.
   3. Global Variables
      a. n_nodes = number of nodes in the graph.
      b. n_edges = number of edges in the graph.
      c. node_list = list of nodes in the graph
   4. Created 5/31/07 by modifying the code from c:\sewell\research\graphs\rgrphgen.c.
      Given the same n_nodes, density, and seed, it should produce the same
      graph as rgrphgen.c.
   5. Modified 3/19/08
      a. Load the information into a bigraph.
      b. Use a character adjacency matrix instead of bit vectors.
      c. Use neighbor_list to store neighbors of nodes.
   6. Modified 12/18/09
      a. Don't load into a bigraph.
      b. Randomly generate the edges, then call build_graph to do the rest.
   7. Modified 1/13/10 to generate node weights.
   8. Modified 3/5/10
      a. Accept a pointer to the graph.
      b. Eliminate weight as an input parameter.
*/
{
   int      i, row, col;
   double   save_seed;

   save_seed = *dseed;

   // Initialize the node names and degrees.

   graph->n_nodes = n;
   //MALLOC(node_list, n_nodes+1, tnode);
   for(i = 1; i <= graph->n_nodes; i++) {
      graph->node_list[i].name = i;
      graph->node_list[i].degree = 0;
   }

   // Initialize the adjacency matrix.

   //MALLOC(adj, n_nodes+1, char *);
   for(row = 1; row <= graph->n_nodes; row++) {
      //MALLOC(adjacent[row], graph->n_nodes+1, char);
      graph->adj[row][row] = 0;
   }

   // Randomly generate the edges.  Compute the degrees.  Create the adjacency matrix.

   graph->n_edges = 0;

   for(row = 1; row <= graph->n_nodes; row++) {
      for(col = row + 1; col <= graph->n_nodes; col++) {
         if(ggubfs(dseed) < density) {
            graph->n_edges = graph->n_edges + 1;
            graph->adj[row][col] = 1;
            graph->adj[col][row] = 1;
         } else {
            graph->adj[row][col] = 0;
            graph->adj[col][row] = 0;
         }
      }
   }

   // Randomly generate the node weights.

   *dseed = save_seed;
   for(i = 1; i <= graph->n_nodes; i++) {
      graph->weight[i] = randomi(11, dseed) - 1;
   }

   build_graph(graph);
}

//_________________________________________________________________________________________________

void testprobs(MWSSgraphpnt graph, MWSSdatapnt data, wstable_parameterspnt parms,
               int n, double density, double* dseed, wstable_infopnt info)
/*
   1. This function randomly generate graphs and solves them.
      an edge being included equals density.
   2. Input Variables
      a. n = number of nodes in the graph.
      b. density = The probability that an edge is created between two nodes.
      c. dseed = seed used for random number generation.
   3. Global Variables
      a. n_nodes = number of nodes in the graph.
      b. n_edges = number of edges in the graph.
      c. node_list = list of nodes in the graph
   4. Written 1/13/10.
   5. Modified 3/5/10
      a. Accept a pointer to the graph.
      b. Eliminate weight as an input parameter.
      c. Accept a pointer to the data structures required for the search.
*/
{
   int      rep;

   for(rep = 1; rep <= 100; rep++) {
      rgrphgen(graph, n, density, dseed);
      initialize_max_wstable(graph, info);
      call_max_wstable(graph, data, parms, info);
      free_reinitialize_graph(graph, data);
   }
   //printf("%10.3f %10.3f\n", info->clique_cover_cpu, info->cpu);
}

//_________________________________________________________________________________________________

double ggubfs(double *dseed)
{
   int      divisor;
   double   product;

      product = 16807.0 * *dseed;
      divisor = product / 2147483647.0;
      *dseed = product - (divisor * 2147483647.0);
      return( *dseed / 2147483648.0 );
}

//_________________________________________________________________________________________________

int randomi(int n, double *dseed)
{
   int      truncate;

   truncate = (n * ggubfs(dseed)) + 1;
   return(truncate);
}

//_________________________________________________________________________________________________

void free_graph(MWSSgraphpnt graph)
{
   int i;
   for(i = 1; i <= graph->n_nodes; i++) {
      //new_memory((void**) &node_list[i].adj_last, n_nodes+1, sizeof(nodepnt *));
      free(graph->node_list[i].adj_last);
   }

   free(graph->edge_list);
   free(graph->adj_last);

}
//_________________________________________________________________________________________________

void free_reinitialize_graph(MWSSgraphpnt graph, MWSSdatapnt data)
{
   int      i, j;
   free_graph(graph);
   for(i = 1; i <= MAX_NODES; i++) {
      graph->node_list[i].active = 0;
      graph->node_list[i].adj = NULL;
      graph->node_list[i].adj2 = NULL;
      if(i <= graph->n_nodes) free(graph->node_list[i].adj_last);
      graph->node_list[i].adj_last = NULL;
      graph->node_list[i].adjacent = NULL;
      graph->node_list[i].adjv = 0;
      graph->node_list[i].degree = 0;
      graph->node_list[i].inverse = 0;
      graph->node_list[i].key = 0;
      graph->node_list[i].name = 0;
      graph->node_list[i].remaining_weight = 0;
      graph->node_list[i].surplus = 0;
      graph->node_list[i].weight = 0;
      graph->weight[i] = 0;

      data->cur_sol[i] = NULL;
      data->n_act[i] = 0;
      data->best_z = 0;
      data->best_sol[i] = NULL;
      data->n_best = 0;
   }

   for(i = 1; i <= MAX_NODES; i++) {
      for(j = 1; j <= MAX_NODES; j++) {
         data->act[i][j] = NULL;
      }
   }
}

//_________________________________________________________________________________________________

int read_dimacs (MWSSgraphpnt graph, char *f, double *weight)
/*
   1. Based on COLORread_dimacs, which was copied from Bill Cook on 12/22/09.
*/
{
      int icount = 0, haveprob = 0;
      int i, j, len, m, n, v, w;
      char buf[256], *p;
      FILE *in = (FILE *) NULL;
      nodepnt  *pn;

      in = fopen (f, "r");
      if (!in) {
         fprintf (stderr, "Unable to open %s for input\n", f);
         goto CLEANUP;
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
               goto CLEANUP;
            }
            haveprob = 1;
            data = strtok(p,delim); /* get 'p' */

            data = strtok(NULL,delim); /* get type */
            if ( strcmp(data,"edge") && strcmp(data,"edges") && strcmp(data,"col") && strcmp(data,"graph") ) {
               fprintf (stderr, "ERROR in Dimacs file -- not an edge file\n");
               goto CLEANUP;
            }
            data = strtok(NULL,delim);
            sscanf (data, "%d", &n);
            graph->n_nodes = n;
            data = strtok(NULL,delim);
            sscanf (data, "%d", &m);
            graph->n_edges = m;

            printf ("Number of Nodes: %d\n", graph->n_nodes);
            printf ("Number of Edges: %d\n", graph->n_edges);
            MALLOC(graph->edge_list, (2*graph->n_edges)+graph->n_nodes, nodepnt);
            MALLOC(graph->adj_last, graph->n_nodes+1, nodepnt*);
            graph->node_list[0].adjacent = NULL;
            for(i = 1; i <= graph->n_nodes; i++) {
               graph->node_list[i].adjacent = graph->adj[i];
            }

            for(i = 0; i <= graph->n_nodes; i++) {
               graph->node_list[i].name = i;
               graph->node_list[i].degree = 0;
               graph->node_list[i].adjv = 0;
               graph->node_list[i].adj2 = NULL;
            }
            graph->node_list[0].name = -1;

            for(i = 0; i <= graph->n_nodes; i++) {
                for(j = 0; j <= graph->n_nodes; j++) {
                  graph->adj[i][j] = 0;
               }
            }

            for(i = 0, pn = graph->edge_list; i < ((2*graph->n_edges)+graph->n_nodes); i++,pn++) {
               *pn = NULL;
            }

            //MALLOC(weight, n_nodes+1, double);
            for (i = 0; i <= graph->n_nodes; i++) weight[i] = 0;
         } else if (p[0] == 'e') {
            if (!haveprob) {
               fprintf (stderr, "ERROR in Dimacs file -- e before p\n");
               goto CLEANUP;
            }
            if (icount >= graph->n_edges) {
               fprintf (stderr, "ERROR in Dimacs file -- to many edges\n");
               goto CLEANUP;
            }
            p++;
            sscanf (p, "%d %d", &v, &w);
            graph->node_list[v].degree++;
            graph->node_list[w].degree++;
            graph->adj[v][w] = 1;
            graph->adj[w][v] = 1;
            icount++;
         } else if (p[0] == 'n') {
            if (!haveprob) {
               fprintf (stderr, "ERROR in Dimacs file -- n before p\n");
               goto CLEANUP;
            }
            p++;
            sscanf (p, "%d %d", &v, &len);
            weight[v] = len;
         }
      }

      pn = graph->edge_list;
      for(i = 1; i <= graph->n_nodes; i++) {
         graph->node_list[i].adj = pn;
         pn += (graph->node_list[i].degree + 1);
      }

      for(i = 1; i < graph->n_nodes; i++) {
         for(j = i+1; j <= graph->n_nodes; j++) {
            if(graph->adj[i][j] == 1) {
               *(graph->node_list[i].adj) = &(graph->node_list[j]);
               if(graph->node_list[i].adj >= pn) {
                  fprintf(stderr, "out of bounds\n");
               }
               graph->node_list[i].adj++;
               *(graph->node_list[j].adj) = &(graph->node_list[i]);
               if(graph->node_list[j].adj >= pn) {
                  fprintf(stderr, "out of bounds\n");
               }
               graph->node_list[j].adj++;
            }
         }
      }

      pn = graph->edge_list;
      for(i = 1; i <= graph->n_nodes; i++) {
         graph->node_list[i].adj = pn;
         graph->adj_last[i] = pn + graph->node_list[i].degree - 1;
         pn += (graph->node_list[i].degree + 1);
      }

CLEANUP:

      if (in) fclose (in);
      return(graph->n_nodes);
 }

