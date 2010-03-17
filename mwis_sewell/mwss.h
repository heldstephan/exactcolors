#ifndef __MWSS_H
#define __MWSS_H
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

#include  <stdio.h>
#include  <stdlib.h>
#include  <time.h>
#include  <assert.h>
#include  <math.h>
#include <string.h>

#define  MALLOC(x,n,type) do {                                         \
         if ((x = (type *) malloc( (n) * sizeof(type))) == NULL) {     \
	         fprintf(stderr,"out of memory\n");                         \
            fprintf(stderr,"x %d type\n",n);                           \
	         exit(1);                                                   \
	 }} while (0)
#define  ABS(i) ((i < 0) ? -(i) : (i) )
#define  MAX(i,j) ((i < j) ? (j) : (i) )
#define  MIN(i,j) ((i < j) ? (i) : (j) )

#define MAX_NODES 1000
#define EPSILON 0.00001

typedef struct node {
   int      name;
   int      degree;
   int      active;
   int      key;                       // sorting key
   int      inverse;                   // index in a sequence
   int      adjv;                      // used by maximal clique to count # of adjacencies
   double   weight;                    // weight of node (used by wstable)
   double   remaining_weight;          // used by weighted_clique_cover
   double   surplus;                   // Weight of v's neighbors minus the weight of v.
   struct   node **adj;
   struct   node **adj2;
   struct   node ***adj_last;
   char     *adjacent;                 // points to the appropriate column of adjacent
} tnode,*nodepnt;

typedef struct wstable_parameters {
   int      clique_cover;              // -k option: clique cover 1 = maximal cliques  2 = maximum cliques
   int      reorder;                   // -r option: sort active nodes by ascending degree
   int      prn_info;                  // Controls level of printed info (def=0)
   double   cpu_limit;                 // cpu limit for search process
} wstable_parameters, *wstable_parameterspnt;

typedef struct wstable_info {
   int      n_subproblems;             // # of subproblems in branch & bound tree
   int      n_sub_depth[MAX_NODES+1];  // # of subproblems at given depth
   int      n_calls;                   // # of times max_wstable is called
   clock_t  start_time;                // starting cpu time
   clock_t  end_time;                  // ending cpu time
   double   cpu;                       // cpu used during search process
   double   clique_cover_cpu;          // cpu used to create clique covers

} wstable_info, *wstable_infopnt;

typedef struct outnode {
   nodepnt  pntv;
   double   w_in_cur_sol;              // Weight of v's neighbors in cur_sol.
} outnode, *outnodepnt;

typedef struct MWSSgraph {
   int      active_flag;                    // used to test active status of a node
   char     adj[MAX_NODES+1][MAX_NODES+1];  // adjacency matrix of the graph
                                            // Note that it does not use first row or column
   int      n_nodes;                // number of nodes in the graph
   int      n_edges;                // number of edges in the graph
   tnode    node_list[MAX_NODES+1]; // vector of nodes of the graph
   nodepnt  *edge_list;             // space for list of edges
   nodepnt  **adj_last;             // points to last active node in adjacency list
   double   weight[MAX_NODES+1];    // weight[v] = weight of node v
} MWSSgraph, *MWSSgraphpnt;

typedef struct MWSSdata {
   nodepnt  cur_sol[MAX_NODES+1];               // Stores the current stable set during search
   nodepnt  act[MAX_NODES+1][MAX_NODES+1];      // Lists of active nodes at each level of search
   int      n_act[MAX_NODES+1];                 // n_active[d] = # of active nodes at depth d
   double   best_z;                             // Weight of best stable set found so far
   nodepnt  best_sol[MAX_NODES+1];              // Stores the best stable set found so far
   int      n_best;                             // # of nodes in best stable set found so far
} MWSSdata, *MWSSdatapnt;

typedef struct osterdata {
   int   oster_cur_sol[MAX_NODES+1];            // Stores the current clique during search
   int   eligible[MAX_NODES+1][MAX_NODES+1];    // Lists of eligible nodes at each level of search
   int   found;                                 // found = 1 if a larger clique/stable set has been found
   int   n_elig[MAX_NODES+1];                   // n_eligible[d] = # of eligible nodes at depth d
   int   oster_best_z;                          // Number of nodes in largest clique/stable set found so far
   int   oster_best_sol[MAX_NODES+1];           // Stores the best clique/stable set found so far
   int   ub[MAX_NODES+1];                       // Used by oster_wclique.
   int   oster_neighbors[MAX_NODES+1];          // Work vector used by maximum_clique
} osterdata, *osterdatapnt;

//extern   double   density;       // -d option: density of graph (randomly generated graph or for printing)
// extern   int      gen;           // -g option: randomly generate problem
// extern   int      lower_bound;   // -l option: user supplied lower bound
// extern   int      n;             // -n option: # of nodes to randomly generate
// extern   double   seed;          // -s option: random seed (def = 3.1567)
// extern   int      test;          // -t option: randomly generate and solve a set of test problems

// Functions in mwss.cpp

void parseargs(int ac, char **av, wstable_parameterspnt parms);
void usage (char *prog);

// Functions in wstable.cpp

void default_parameters(wstable_parameterspnt parms);
void initialize_max_wstable(MWSSgraphpnt graph, wstable_infopnt info);
void free_graph(MWSSgraphpnt graph);

void call_max_wstable(MWSSgraphpnt graph, MWSSdatapnt data,
                      wstable_parameterspnt parameters,
                      wstable_infopnt info,
		      double goal);
int max_wstable(MWSSgraphpnt graph, MWSSdatapnt data, nodepnt *best_stable, int *n_best_stable, double *z_best, wstable_infopnt info,
                wstable_parameterspnt parameters, nodepnt *list, int n_list, double lower_bound,
                double goal);
void wstable(MWSSgraphpnt graph, MWSSdatapnt data, osterdatapnt oster_data, int adj_last_offset, double cur_z, int depth, double goal, outnode *out, int n_out,
           wstable_infopnt info, wstable_parameterspnt parameters);
void determine_branch_nodes(MWSSgraphpnt graph, MWSSdatapnt data, osterdatapnt oster_data, nodepnt *active, int adj_last_offset, double cur_z, int n_active,
                       int *n_branch_nodes, outnode *out, int n_out, wstable_infopnt info, wstable_parameterspnt parameters);
void stable_sub_problem(MWSSdatapnt data, int adj_last_offset, int *branch, nodepnt branch_node, double cur_z,
                        int depth, int j, int *n_sub_active);
void weighted_clique_cover(MWSSgraphpnt graph, osterdatapnt oster_data, nodepnt *active, nodepnt *active2, int adj_last_offset, nodepnt *branch_nodes,
              int n_active, int *n_active2, double available_weight, int *n_branch_nodes, wstable_infopnt info, wstable_parameterspnt parameters);
void maximal_wclique(MWSSgraphpnt graph, int adj_last_offset, nodepnt pntv);
void maximum_clique(MWSSgraphpnt graph, osterdatapnt oster_data, int adj_last_offset, nodepnt pntv);
void oster_clique(MWSSgraphpnt graph, osterdatapnt oster_data, int cur_z, int depth, char clique_or_stable);
int oster_maximum_clique(MWSSgraphpnt graph, osterdatapnt oster_data, int list[], int n_list, char clique_or_stable);
int check_cliq(MWSSgraphpnt graph, osterdatapnt oster_data, char clique_or_stable, int n_best);
void branch_surplus(nodepnt *active, nodepnt *active2, nodepnt *branch_nodes, int n_active,
                    int *n_active2, int *n_branch_nodes);
void branch_neighbors(nodepnt *active, nodepnt *active2, /* int adj_last_offset, */ nodepnt *branch_nodes,
                      int n_active, int *n_active2, int *n_branch_nodes);
void reorder_nodes(MWSSgraphpnt graph, int adj_last_offset, nodepnt *list, int n_list);
void reorder_edges(MWSSgraphpnt graph, int adj_last_offset, nodepnt *list, int n_list);
void branch_out(nodepnt *active, nodepnt *active2, nodepnt *branch_nodes, int n_active,
                int *n_active2, int *n_branch_nodes, outnode *out, int n_out);
int reorder_out(MWSSgraphpnt graph, int adj_last_offset/* , int depth */, int n_out, outnode *out);
void new_out(int adj_last_offset, nodepnt branch_node, outnode *out, int n_out, outnode *sub_out, int *n_sub_out);
void ascending_distrib_sort(nodepnt *list, int m, int n);
double greedy_wstable(MWSSgraphpnt graph, nodepnt *list, int n_list, nodepnt *stable_set, int *n_stable_set);
void build_graph(MWSSgraphpnt graph);
int check_graph(MWSSgraphpnt graph);
int check_subgraph(MWSSgraphpnt graph, int adj_last_offset, nodepnt *list, int n_list);
void prn_graph(MWSSgraphpnt graph);
void prn_subgraph(int adj_last_offset, nodepnt *list, int n_list);
void prn_nodes(nodepnt *list, int n);
void prn_node_weights(nodepnt *list, int n);
void prn_stable(nodepnt *active, double best_z, double cur_z, nodepnt *branch_nodes, int depth,
                int n_active, int n_branch_nodes, wstable_infopnt info, int prn_level);
void rgrphgen(MWSSgraphpnt graph, int n, double density, double* dseed);
void testprobs(MWSSgraphpnt graph, MWSSdatapnt data, wstable_parameterspnt parms,
               int n, double density, double* dseed, wstable_infopnt info);
double ggubfs(double *dseed);
int randomi(int n, double *dseed);
void free_reinitialize_graph(MWSSgraphpnt graph, MWSSdatapnt data);
int read_dimacs (MWSSgraphpnt graph, char *f, double *weight);

#endif
