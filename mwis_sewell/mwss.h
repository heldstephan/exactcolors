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

#define  MWIS_MALLOC(x,n,type) do {                                      \
      x = (type *) malloc( (n) * sizeof(type));                          \
   } while (0)

#define  MWIS_IFFREE(x,type) do {                                         \
      if (x) free(x);                                                     \
      x = (type*) NULL;                                                   \
   } while (0)

#define MWIScheck_NULL(item,msg) {                                        \
    if ((!item)) {                                                         \
       fflush(stdout);                                                     \
       fprintf (stderr, "%s at %s, line %d\n", (msg),__FILE__,__LINE__);   \
       rval = 1;                                                           \
       goto CLEANUP;                                                       \
    }                                                                      \
}
#define MWIScheck_rval(rval,msg) {                                        \
    if ((rval)) {                                                          \
       fflush(stdout);                                                     \
       fprintf (stderr, "%s at %s, line %d\n", (msg),__FILE__,__LINE__);   \
       goto CLEANUP;                                                       \
    }                                                                      \
}



#define  ABS(i) ((i < 0) ? -(i) : (i) )
#define  MAX(i,j) ((i < j) ? (j) : (i) )
#define  MIN(i,j) ((i < j) ? (i) : (j) )


typedef int MWISNW;

#define MWISNW_EPSILON 0
#define MWISNW_MAX     INT_MAX

typedef struct node {
   int      name;
   int      degree;
   int      active;
   int      key;                       // sorting key
   int      inverse;                   // index in a sequence
   int      adjv;                      // used by maximal clique to count # of adjacencies
   MWISNW   weight;                    // weight of node (used by wstable)
   MWISNW   remaining_weight;          // used by weighted_clique_cover
   MWISNW   surplus;                   // Weight of v's neighbors minus the weight of v.
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
   int*     n_sub_depth;               // # of subproblems at given depth
   int      n_calls;                   // # of times max_wstable is called
   clock_t  start_time;                // starting cpu time
   clock_t  end_time;                  // ending cpu time
   double   cpu;                       // cpu used during search process
   double   clique_cover_cpu;          // cpu used to create clique covers

} wstable_info, *wstable_infopnt;

typedef struct outnode {
   nodepnt  pntv;
   MWISNW   w_in_cur_sol;              // Weight of v's neighbors in cur_sol.
} outnode, *outnodepnt;

typedef struct MWSSgraph {
   int      active_flag;                    // used to test active status of a node
   char**   adj;                            // adjacency matrix of the graph
                                            // Note that it does not use first row or column
   int      n_nodes;                // number of nodes in the graph
   int      n_edges;                // number of edges in the graph
   tnode*   node_list;              // vector of nodes of the graph
   nodepnt  *edge_list;             // space for list of edges
   nodepnt  **adj_last;             // points to last active node in adjacency list
   MWISNW*   weight;                // weight[v] = weight of node v
} MWSSgraph, *MWSSgraphpnt;

typedef struct MWSSdata {
   nodepnt*  cur_sol;                           // Stores the current stable set during search
   nodepnt** act;                               // Lists of active nodes at each level of search
   int*      n_act;                             // n_active[d] = # of active nodes at depth d
   MWISNW   best_z;                             // Weight of best stable set found so far
   nodepnt* best_sol;                           // Stores the best stable set found so far
   int      n_best;                             // # of nodes in best stable set found so far
   int      n_nodes;                            // number of nodes in graph.
} MWSSdata, *MWSSdatapnt;

typedef struct osterdata {
   int*  oster_cur_sol;            // Stores the current clique during search
   int** eligible;                 // Lists of eligible nodes at each level of search
   int   found;                                 // found = 1 if a larger clique/stable set has been found
   int*  n_elig;                   // n_eligible[d] = # of eligible nodes at depth d
   int   oster_best_z;                          // Number of nodes in largest clique/stable set found so far
   int*  oster_best_sol;           // Stores the best clique/stable set found so far
   int*  ub;                       // Used by oster_wclique.
   int*  oster_neighbors;          // Work vector used by maximum_clique
   int   n_nodes;
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

int initialize_max_wstable(MWSSgraphpnt graph, wstable_infopnt info);
void free_max_wstable(MWSSgraphpnt graph, MWSSdatapnt data, wstable_infopnt info);

int call_max_wstable(MWSSgraphpnt graph, MWSSdatapnt data,
                     wstable_parameterspnt parameters,
                     wstable_infopnt info,
                     MWISNW goal);
int max_wstable(MWSSgraphpnt graph, MWSSdatapnt data, nodepnt *best_stable, int *n_best_stable, MWISNW *z_best, wstable_infopnt info,
                wstable_parameterspnt parameters, nodepnt *list, int n_list, MWISNW lower_bound,
                MWISNW goal, int* status);
int wstable(MWSSgraphpnt graph, MWSSdatapnt data, osterdatapnt oster_data, int adj_last_offset, MWISNW cur_z, int depth, MWISNW goal, outnode *out, int n_out,
            wstable_infopnt info, wstable_parameterspnt parameters);
int determine_branch_nodes(MWSSgraphpnt graph, MWSSdatapnt data, osterdatapnt oster_data, nodepnt *active, int adj_last_offset, MWISNW cur_z, int n_active,
                           int *n_branch_nodes, outnode *out, int n_out, wstable_infopnt info, wstable_parameterspnt parameters);
void stable_sub_problem(MWSSdatapnt data, int adj_last_offset, int *branch, nodepnt branch_node, MWISNW cur_z,
                        int depth, int j, int *n_sub_active);
void weighted_clique_cover(MWSSgraphpnt graph, osterdatapnt oster_data, nodepnt *active, nodepnt *active2, int adj_last_offset, nodepnt *branch_nodes,
              int n_active, int *n_active2, MWISNW available_weight, int *n_branch_nodes, wstable_infopnt info, wstable_parameterspnt parameters);
int maximal_wclique(MWSSgraphpnt graph, int adj_last_offset, nodepnt pntv);
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
int ascending_distrib_sort(nodepnt *list, int m, int n);
int greedy_wstable(MWSSgraphpnt graph, nodepnt *list, int n_list, nodepnt *stable_set, int *n_stable_set, 
                   MWISNW* greedy_weight);
void reset_pointers(MWSSgraphpnt graph, MWSSdatapnt data, wstable_infopnt info);
void reset_osterdata_pointers(osterdatapnt odata);

int allocate_data(MWSSdatapnt data, int n_nodes);
int allocate_graph(MWSSgraphpnt graph, int n_nodes);
int allocate_osterdata(osterdatapnt oster_data, int n_nodes);

void free_graph(MWSSgraphpnt graph);
void free_data(MWSSdatapnt data);
void free_wstable_info(wstable_infopnt info);
void free_osterdata(osterdatapnt oster_data);


int build_graph(MWSSgraphpnt graph);
int check_graph(MWSSgraphpnt graph);
int check_subgraph(MWSSgraphpnt graph, int adj_last_offset, nodepnt *list, int n_list);
void prn_graph(MWSSgraphpnt graph);
void prn_subgraph(int adj_last_offset, nodepnt *list, int n_list);
void prn_nodes(nodepnt *list, int n);
void prn_node_weights(nodepnt *list, int n);
void prn_stable(nodepnt *active, MWISNW best_z, MWISNW cur_z, nodepnt *branch_nodes, int depth,
                int n_active, int n_branch_nodes, wstable_infopnt info, int prn_level);
int rgrphgen(MWSSgraphpnt graph, int n, double density, double* dseed);
int testprobs(MWSSgraphpnt graph, MWSSdatapnt data, wstable_parameterspnt parms,
              int n, double density, double* dseed, wstable_infopnt info);
double ggubfs(double *dseed);
int randomi(int n, double *dseed);
void free_reinitialize_graph(MWSSgraphpnt graph, MWSSdatapnt data);
int read_dimacs (MWSSgraphpnt graph, char *f);

#endif
