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
 * dsatur.c arises from:
 * Janina Lea Vogl: Exakte Algorithmen für Graphenfärbungen, bachelor's thesis, University of Bonn, 2018.
 *
 * Die Funtion DSATUR ist die Hauptfunktion dieser Datei, die den modifizierten DSATUR-Algorithmus implementiert,
 * der die fraktionale chromatische Zahl verwendet, um die lower_bound zu verbessern.
 * DSATUR_recursion ist die zugehörige Rekursionsfunktion
 * Alle Funktionen mit DSATUR_ am Anfang sind selbstgeschriebene/abgewandelte Hilfsfunktionen für DSATUR oder DSATUR_recursion
 * print_COLORproblem_edges, print_COLORadjgraph, copy_colors sind auch selbstgeschrieben 
 * restliche Funktionen sind kopierte Funktionen, die nicht in header-Dateien gelistet sind, deshalb in der Datei stehen
 * Kompilieren mit dem dsaturMakefile
 */

#include "color.h"
#include "graph.h"
#include "color_defs.h"
#include "color_private.h"
#include <assert.h>

unsigned long nodes=0;

typedef struct DSATURREC
{
    int color_count;	//Anzahl der Farben
    int global_lower_bound;	//globale untere Schranke LB
    int uncolored_nodes_number;	//Anzahl ungefärbter Knoten
    int optimal;	//0, wenn LB < UB und 1, wenn LB=UB => bricht Rekursion ab
    int *DSAT_value;	//enthält Anzahl verschiedener Farben in Nachbarschaft für jeden Knoten
    int *id_in_child_cd;	//Index: Id in G und Eintrag: Id in child_cd
} DSATURREC;

//Hilfsstruct für die Sortierung der Knoten in DSATUR_sort_vertices_by_degree
typedef struct sort
{
    int id;        //id des Knoten im Ursprungsgraphen
    int in_clique; //1 wenn Knoten in der Clique, sonst 0
    int degree;    //Grad des Knotens
} sort;

static int DSATUR_create_differ(colordata* parent_cd, int v1, int v2);
int DSATURREC_init(DSATURREC *dsat, int lowerbound);
void print_COLORproblem_edges(colordata *p);
void print_COLORadjgraph(COLORadjgraph *G);
int DSATUR_simplify(COLORadjgraph *G, COLORadjgraph *H, int **node_mapping, int global_lower_bound);
int copy_colors(COLORadjgraph *G, int color_count, COLORset **sets, int *count);
void DSATUR_update_dsat_value_for_neighbors(DSATURREC *dsat, int colored_node, COLORadjgraph *G, int *feasible_colors, int *updated_neighbors);
int DSATUR_create_branch(colordata *cd, COLORproblem *problem, COLORadjgraph *G, int v1, int v2);
static int DSATUR_comp_by_degree(const void *v1, const void *v2);
int DSATUR_sort_vertices_by_degree(COLORadjgraph *G, int **node_mapping, int n, int *clique, int clique_number);
int DSATUR_MaxClique_recursion(int **clique, int *max, int *found, int *clique_sizes_for_nodes, int *marked_nodes, COLORadjgraph *G, int **set, int *set_size, int size);
int DSATUR_MaxClique(int **clique, COLORadjgraph *G, int *lower_bound);
int DSATUR_choose_and_compute_lower_bound(DSATURREC *dsat,int **clique, COLORadjgraph *G, COLORproblem *p);
int DSATUR_new_COLORproblem(COLORproblem *new_problem, COLORproblem *problem, COLORadjgraph *G, int color_count);
int DSATUR(COLORproblem *problem, int *ncolors, COLORset **colorclasses);
int DSATUR_recursion(colordata *child_cd, COLORproblem *problem, COLORadjgraph *G, DSATURREC *dsat, int *feasible_colors, int* aufrufanzahl, double *lbzeit, int *colored_node_updated_neighbors);

//Funktionen aus bereits vorhandenen Dateien, die nicht durch header-Dateien erreichbar sind, mit einigen weiteren Kommentaren
static int comp_node_ids(const void* v1, const void* v2)
{
   int id1 = * (const int*) v1;
   int id2 = * (const int*) v2;

   return id1 - id2;
}

//färbt Knoten in minimaler zulässiger Farbe
static void color_node (COLORadjgraph *G, int n)
{
    int i, color = 0;
    COLORadjnode *p = &G->nodelist[n];
    int failed;
    do {    //geht solange Farben durch, bis die Farbe nicht mehr in der Nachbarschaft vorkommt
     failed = 0;
       for (i = 0; !failed && i < p->degree; i++) { //kommt diese Farbe in der Nachbarschaft vor?
          if (G->nodelist[p->adj[i]].color == color) {
             ++color;
             failed = 1;
          }
       }
    } while (failed);
    p->color = color;
}

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
    fprintf (stderr, "   -s int Branching strategy: 0 = none, 1 = minimum lower bound (default),"
             " 2 = DFS, 3 = hybrid (2 followed by 1).\n");
    fprintf (stderr, "   -R int rounding style: 0 = neighbor (default), 1 = uniform, 2 = none\n");
    fprintf (stderr, "   -l dbl cpu time limit for branching.\n");

}

//setzt in parms, entsprechend der Kommandozeilenparameter, die Parameter
static int parseargs(int ac, char **av, COLORparms *parms)
{
    int c;
    int rval = 0;
    int debug = COLORdbg_lvl();

    while ((c = getopt(ac, av, "admpo:r:w:c:u:b:l:s:R:")) != EOF)
    {
        switch (c)
        {
        case 'd':
            /* each -d increases the verbosity by one.*/
            ++(debug);
            COLORset_dbg_lvl(debug);
            break;
        case 'o':
            rval = COLORparms_set_outfile(parms, optarg);
            COLORcheck_rval(rval, "Failed in COLORparms_set_outfile");
            break;
        case 'm':
            COLORparms_set_write_mwis(parms, 1);
            break;
        case 'r':
            rval = COLORparms_set_cclasses_infile(parms, optarg);
            COLORcheck_rval(rval, "Failed in COLORparms_set_cclasses_infile");
            break;
        case 'w':
            rval = COLORparms_set_cclasses_outfile(parms, optarg);
            COLORcheck_rval(rval, "Failed in COLORparms_set_cclasses_outfile");
            break;
        case 'c':
            rval = COLORparms_set_color_infile(parms, optarg);
            COLORcheck_rval(rval, "Failed in COLORparms_set_color_infile");
            break;
        case 'u':
            rval = COLORparms_set_initial_upper_bound(parms, atoi(optarg));
            COLORcheck_rval(rval, "Failed in COLORparms_set_initial_upper_bound");
            break;
        case 'a':
            parms->upper_bounds_only = 1;
            parms->branching_strategy = COLOR_hybrid_strategy;
            break;
        case 'p':
            rval = COLORparms_set_parallel(parms, 1);
            COLORcheck_rval(rval, "Failed in COLORparms_set_initial_upper_bound");
            break;
        case 'b':
            rval = COLORparms_set_backupdir(parms, optarg);
            COLORcheck_rval(rval, "Failed in COLORparms_set_backupdir");
            break;
        case 'l':
            rval = COLORparms_set_branching_cpu_limit(parms, atof(optarg));
            COLORcheck_rval(rval, "Failed in COLORparms_set_branching_cpu_limit");
            break;
        case 's':
            rval = COLORparms_set_branching_strategy(parms, atoi(optarg));
            COLORcheck_rval(rval, "Failed in COLORparms_set_branching_strategy");
            break;
        case 'R':
            rval = COLORparms_set_rounding_strategy(parms, atoi(optarg));
            COLORcheck_rval(rval, "Failed in COLORparms_set_rounding_strategy");
            break;
        default:
            usage(av[0]);
            rval = 1;
            goto CLEANUP;
        }
    }

    if (ac <= optind)
    {
        rval = 1;
        goto CLEANUP;
    }
    else
    {
        rval = COLORparms_set_edgefile(parms, av[optind++]);
        COLORcheck_rval(rval, "Failed in COLORparms_set_edgefile");
    }

CLEANUP:

    if (rval)
        usage(av[0]);
    return (rval);
}

//gibt Problemnamen aus und speichert ihn in pname
static int get_problem_name(char *pname, const char *efname)
{
    int rval = 0;
    int len = 0;
    const char *fname = strrchr(efname, '/');
    char *lastdot = strrchr(efname, '.');
    if (!fname)
    {
        /* no slashes in efname.*/
        fname = efname;
    }
    else
    {
        fname++;
    }

    if (lastdot)
    {
        len = lastdot - fname + 1;
    }
    else
    {
        len = strlen(fname);
    }

    if (snprintf(pname, len, "%s", fname) < 0)
    {
        rval = 1;
    }
    printf("Extracted problem name %s\n", pname);

    return rval;
}


static void COLORset_SWAP(COLORset *c1,COLORset *c2,COLORset *t)
{
   if (c1 != c2) {
      memcpy(t,c2,sizeof(COLORset));
      memcpy(c2,c1,sizeof(COLORset));
      memcpy(c1,t,sizeof(COLORset));
   }
}

static int COLORset_less(COLORset *c1,COLORset *c2)
{
   int i;
   if (c1->count != c2->count) {
      return c1->count < c2->count;
   }
   for (i = 0; i < c1->count;++i) {
      if (c1->members[i] != c2->members[i]) {
         return c1->members[i] < c2->members[i];
      }
   }
   return 0;
}

static void COLORset_unify (COLORset *cclasses, int* new_ccount, int ccount)
{
   int i;
   COLORset temp;
   COLORset_quicksort(cclasses,ccount);

   *new_ccount = 0;
   i = 0;
   if (! ccount) return;

   /* Find first non-empty set */
   while(!cclasses[i].count) {
      COLORfree_set(&(cclasses[i++]));
   }

   for (; i < ccount; ++i) {
      if (*new_ccount == 0 || COLORset_less(&(cclasses[*new_ccount -1]),&(cclasses[i]))) {
         (*new_ccount)++;
         if(*new_ccount  < i + 1) {
            COLORset_SWAP(&(cclasses[*new_ccount -1]),& (cclasses[i]),&temp);
         }
      } else {
         COLORfree_set(&(cclasses[i]));
      }
   }
}


static int new_eindex(int v, int v1, int v2)
{
   if (v == v2) {return v1;}
   if (v > v2)  {return v - 1;}

   return v;
}

static int mark_neighborhood(int* neighbor_marker,
                             int ncount, int ecount, int elist[],
                             int v)

{
   int rval = 0;
   int i;
   COLORadjgraph G;

   COLORadjgraph_init(&G);
   rval = COLORadjgraph_build(&G,ncount,ecount, elist);
   COLORcheck_rval(rval,"COLORadjgraph_build");

   for(i = 0; i < ncount; ++i) { neighbor_marker[i] = 0;}

   for(i = 0; i < G.nodelist[v].degree; ++i) {
      int v_i = G.nodelist[v].adj[i];
      neighbor_marker[v_i] = 1;
   }
   COLORadjgraph_free(&G);

 CLEANUP:
   return rval;
}

static int contract_elist(int *elist[],
                          int ncount, int* ecount,
                          int v1, int v2)
{
   int rval = 0;
   int i;
   COLORadjgraph G;

   for(i = 0; i < *ecount; ++i) {
      (*elist)[2*i]   = new_eindex((*elist)[2*i],v1,v2);
      (*elist)[2*i+1] = new_eindex((*elist)[2*i+1],v1,v2);
   }

   COLORadjgraph_init(&G);
   rval = COLORadjgraph_build(&G, ncount,*ecount, *elist);
   COLORcheck_rval(rval,"COLORadjgraph_build");

   rval =  COLORadjgraph_simplify(&G);
   COLORcheck_rval(rval,"COLORadjgraph_simplify");

   rval = COLORadjgraph_extract_edgelist(ecount, elist,&G);

 CLEANUP:
   COLORadjgraph_free(&G);

   return rval;
}

static int create_contracted_graph(colordata* cd,
                                   int ncount, int ecount, int elist[],
                                   int v1, int v2)
{
   int rval = 0;
   int i;


   /* Create contracted graph */
   cd->ncount = ncount - 1;
   cd->ecount = ecount;

   cd->elist  = (int*) COLOR_SAFE_MALLOC(2 * cd->ecount,int);
   COLORcheck_NULL(cd->elist,"Failed to allocate cd->elist");

   cd->orig_node_ids = (int*) COLOR_SAFE_MALLOC(cd->ncount,int);
   COLORcheck_NULL(cd->orig_node_ids,"Failed to allocate cd->orig_node_ids");

   memcpy(cd->elist,elist, 2 * ecount * sizeof(int));

   rval = contract_elist(&cd->elist,cd->ncount,&cd->ecount,v1,v2);

   for (i = 0; i < cd->ncount; ++i) {
      int j = (i < v2 ) ?  i : i + 1;
      cd->orig_node_ids[i] = cd->parent->orig_node_ids[j];
   }

 CLEANUP:

   return rval;
}

static int prune_duplicated_sets (colordata* cd)
{
   int rval = 0;
   int i,j;
   COLORset_unify (cd->cclasses, &(cd->ccount), cd->ccount);
   for (i = 0 ; i < cd->ccount; ++i ) {
      if (COLORdbg_lvl() > 1) {
            printf("TRANSSORT SET ");
            for (j = 0; j < cd->cclasses[i].count; ++j) {
               printf(" %d",cd->cclasses[i].members[j]);
            }
            printf("\n");
      }
      rval = COLORcheck_set(&(cd->cclasses[i]),cd->ncount,cd->ecount,cd->elist);
      COLORcheck_rval(rval, "Illegal colorset created in create_same\n!");
   }
 CLEANUP:
   return rval;
}

static int transfer_same_cclasses(colordata* cd,
                                  const int* v2_neighbor_marker,
                                  const COLORset* parent_cclasses,
                                  int   parent_ccount,
                                  int   v1,
                                  int   v2)
{
   int rval = 0;
   int i;
   int v1_covered = 0;
   /* Transfer independent sets: */
   cd->gallocated = cd->ccount   =  parent_ccount + 1;
   cd->cclasses = (COLORset*) COLOR_SAFE_MALLOC(cd->gallocated,COLORset);
   for (i = 0; i < parent_ccount; ++i) {
      int j;
      int add_v1 = 1;

      COLORinit_set(cd->cclasses + i);

      cd->cclasses[i].members = (int*) COLOR_SAFE_MALLOC(parent_cclasses[i].count,int);
      cd->cclasses[i].count = 0;
      for (j = 0; j < parent_cclasses[i].count; ++j) {
         if (v2_neighbor_marker[parent_cclasses[i].members[j]] == 1) {
            add_v1 = 0;
            j = parent_cclasses[i].count;/*break*/
         }
      }
      for (j = 0; j < parent_cclasses[i].count; ++j) {
         if (parent_cclasses[i].members[j] == v1) {
            if (add_v1) {v1_covered = 1;}
            else {continue;}
         }
         if (parent_cclasses[i].members[j] < v2) {
            cd->cclasses[i].members[(cd->cclasses[i].count)++] =
               parent_cclasses[i].members[j];
         }
         else if (parent_cclasses[i].members[j] >  v2) {
            cd->cclasses[i].members[(cd->cclasses[i].count)++] =
               parent_cclasses[i].members[j] - 1;
         }
         /* else 'parent_cclasses[i].members[j] == v2' and we skip it*/
      }
      if(COLORdbg_lvl() > 1) {
         printf("PARENT SET ");
         for (j = 0; j < parent_cclasses[i].count; ++j) {
            printf(" %d",parent_cclasses[i].members[j]);
         }
         printf("\n");
         printf("TRANS SET ");
         for (j = 0; j < cd->cclasses[i].count; ++j) {
            printf(" %d",cd->cclasses[i].members[j]);
         }
         printf("\n");
      }
   }
   if (!v1_covered) {
      /* Create the singular set v1 as last set, because we might not add
         v1 to any set: */
      printf("Adding extra set %d\n", v1);
      cd->cclasses[parent_ccount].count  = 1;
      cd->cclasses[parent_ccount].members = (int*) COLOR_SAFE_MALLOC(1,int);
      cd->cclasses[parent_ccount].members[0] = v1;
   } else {
      cd->cclasses[parent_ccount].count  = 0;
      cd->cclasses[parent_ccount].members = (int*) 0;
      cd->ccount--;
   }
   rval = prune_duplicated_sets(cd);
   COLORcheck_rval(rval,"Failed in prune_duplicated_sets");

   /* END Transfer independent sets: */

 CLEANUP:
   return rval;
}

static int create_same (colordata* parent_cd,
                        int v1, int v2)
{
   int rval = 0;
   colordata*    cd = (colordata*) NULL;
   int* v2_neighbor_marker = (int*) COLOR_SAFE_MALLOC(parent_cd->ncount,int);
   COLORcheck_NULL(v2_neighbor_marker,"Failed to allocate v2_neighbor_marker");

   cd = (colordata*) COLOR_SAFE_MALLOC(1,colordata);
   COLORcheck_NULL(cd,"Failed to allocate cd");

   init_colordata(cd);
   cd->depth = parent_cd->depth + 1;
   parent_cd->nsame         = 1;
   parent_cd->same_children = cd;

   cd->v1 = v1;
   cd->v2 = v2;

   cd->upper_bound = parent_cd->upper_bound;
   cd->lower_bound = parent_cd->lower_bound;
   cd->dbl_safe_lower_bound = parent_cd->dbl_safe_lower_bound;

   cd->parent = parent_cd;
   cd->debugcolors = parent_cd->debugcolors;
   cd->ndebugcolors = parent_cd->ndebugcolors;

   rval = mark_neighborhood(v2_neighbor_marker,
                            parent_cd->ncount,parent_cd->ecount,parent_cd->elist,
                            v2);
   COLORcheck_rval(rval,"Failed in mark_neighborhood");

   rval = create_contracted_graph (cd,
                                   parent_cd->ncount,
                                   parent_cd->ecount,parent_cd->elist,
                                   v1,v2);
   COLORcheck_rval(rval,"Failed in create_contracted_graph_and_intersect");


   /* END create contracted graph */

   if (COLORdbg_lvl() > 1) {
      printf("create_same created following graph:\n");
      COLORgraph_print(cd->ecount,cd->elist);
   }

   rval = transfer_same_cclasses(cd,
                                 v2_neighbor_marker,
                                 parent_cd->cclasses,
                                  parent_cd->ccount,
                                 v1,v2);
   COLORcheck_rval(rval,"Failed in transfer_same_cclasses");
 CLEANUP:
   if (rval) {
      if (cd) {
         free_colordata(cd);
         free(cd);
      }
      parent_cd->same_children = (colordata*) NULL;
   }
   if (v2_neighbor_marker) free(v2_neighbor_marker);
   return rval;
}



//Modifizierte create_differ Funktion, damit auch Sonderfall v1>=v2 funktioniert
static int DSATUR_create_differ(colordata* parent_cd, int v1, int v2)
{
   int rval = 0;
   int i;
   COLORadjgraph G;
   int           v2_covered = 0;
   colordata*    cd = (colordata*) COLOR_SAFE_MALLOC(1,colordata);
   COLORcheck_NULL(cd,"Failed to allocate cd");

   init_colordata(cd);

   cd->depth = parent_cd->depth + 1;
   parent_cd->ndiff         = 1;
   parent_cd->diff_children = cd;

   cd->v1 = v1;
   cd->v2 = v2;

   cd->upper_bound = parent_cd->upper_bound;
   cd->lower_bound = parent_cd->lower_bound;
   cd->dbl_safe_lower_bound = parent_cd->dbl_safe_lower_bound;

   /* Create  graph with extra edge (v1,v2) */
   cd->ncount = parent_cd->ncount;
   cd->ecount = parent_cd->ecount;
   if(v1 != v2)
   {
     cd->ecount++;
   }
   cd->elist  = (int*) COLOR_SAFE_MALLOC(2 * cd->ecount,int);
   COLORcheck_NULL(cd->elist,"Failed to allocate cd->elist");

   memcpy(cd->elist,parent_cd->elist, 2 * parent_cd->ecount * sizeof(int));
   if(v1 != v2){
      cd->elist[ 2 * (cd->ecount - 1)] = v1;
      cd->elist[ 2 * (cd->ecount - 1) + 1] = v2;
    
      rval = COLORadjgraph_build(&G, cd->ncount,cd->ecount,cd->elist);
      COLORcheck_rval(rval,"COLORadjgraph_build failed");

      COLORadjgraph_sort_adjlists_by_id(&G);

      COLORadjgraph_extract_edgelist(&cd->ecount, &cd->elist,&G);
      COLORcheck_rval(rval,"COLORadjgraph_extract_edgelist failed");
   
      /* END: Create  graph with extra edge (v1,v2) */

      if (COLORdbg_lvl() > 1) {
	  printf("DSATUR_create_differ created following graph:\n");
	  COLORgraph_print(cd->ecount,cd->elist);
      }
   }
   
   cd->orig_node_ids = (int*) COLOR_SAFE_MALLOC(cd->ncount,int);
   COLORcheck_NULL(cd->orig_node_ids,"Failed to allocate cd->orig_node_ids");
   for (i = 0; i < cd->ncount; ++i) {
      cd->orig_node_ids[i] = parent_cd->orig_node_ids[i];
   }
   cd->parent = parent_cd;
   cd->debugcolors = parent_cd->debugcolors;
   cd->ndebugcolors = parent_cd->ndebugcolors;
   
   /* Transfer independent sets by removing v2 if both v1 and v2 are currently contained: */
   cd->gallocated = cd->ccount   =  parent_cd->ccount + 1;
   cd->cclasses = (COLORset*) COLOR_SAFE_MALLOC(cd->gallocated,COLORset);

   COLORinit_set(cd->cclasses + parent_cd->ccount);

   for (i = 0; i < parent_cd->ccount; ++i) {
      int j;
      int v1_found = 0;

      COLORinit_set(cd->cclasses + i);

      cd->cclasses[i].members = (int*) COLOR_SAFE_MALLOC(parent_cd->cclasses[i].count,int);
      COLORcheck_NULL(cd->cclasses[i].members,"Failed to allocate cd->cclasses[i].members");
      cd->cclasses[i].count = 0;

      if(v1 < v2){
	  for (j = 0; j < parent_cd->cclasses[i].count; ++j) {
	    int current_elm = parent_cd->cclasses[i].members[j];
	    if (current_elm ==  v1) {
		v1_found = 1;
	    }
	    if (current_elm ==  v2) {
		if (v1_found) {
		  continue;
		} else {
		  v2_covered = 1;
		}
	    }
	    cd->cclasses[i].members[cd->cclasses[i].count] = current_elm;
	    (cd->cclasses[i].count)++;
	  }
      }  
      else if(v1 > v2){
	  for(j = 1; j < parent_cd->cclasses[i].count;++j){
	      if(parent_cd->cclasses[i].members[j] ==  v1) {
		  v1_found = 1;
		  break;
	      }
	  }
	  for (j = 0; j < parent_cd->cclasses[i].count; ++j) {
	    int current_elm = parent_cd->cclasses[i].members[j];
	    if (current_elm ==  v2) {
		if (v1_found) {
		  continue;
		} else {
		  v2_covered = 1;
		}
	    }
	    cd->cclasses[i].members[cd->cclasses[i].count] = current_elm;
	    (cd->cclasses[i].count)++;
	  }
      }
      else{
	  for (j = 0; j < parent_cd->cclasses[i].count; ++j) {
	    int current_elm = parent_cd->cclasses[i].members[j];
	    if (current_elm ==  v2) {
		v2_covered = 1;
	    }
	    cd->cclasses[i].members[cd->cclasses[i].count] = current_elm;
	    (cd->cclasses[i].count)++;
	  }
      }
      if(COLORdbg_lvl() > 1 ) {
         printf("PARENT SET ");
         for (j = 0; j < parent_cd->cclasses[i].count; ++j) {
            printf(" %d",parent_cd->cclasses[i].members[j]);
         }
         printf("\n");

         printf("TRANS SET ");
         for (j = 0; j < cd->cclasses[i].count; ++j) {
            printf(" %d",cd->cclasses[i].members[j]);
         }
         printf("\n");
      }
      rval = COLORcheck_set(&(cd->cclasses[i]),cd->ncount,cd->ecount,cd->elist);
      COLORcheck_rval(rval, "Illegal colorset created in DSATUR_create_differ\n!");

   }
   if (!v2_covered) {
      /* Create the singular set v2 as last set, because we might not add
         v2 to any set: */
      cd->cclasses[parent_cd->ccount].count  = 1;
      cd->cclasses[parent_cd->ccount].members = (int*) COLOR_SAFE_MALLOC(1,int);
      cd->cclasses[parent_cd->ccount].members[0] = v2;
   } else {
      cd->cclasses[parent_cd->ccount].count   = 0;
      cd->cclasses[parent_cd->ccount].members = (int*) 0;
      cd->ccount--;
   }
   /* END Transfer independent sets: */


   rval = prune_duplicated_sets(cd);
   COLORcheck_rval(rval,"Failed in prune_duplicated_sets");


 CLEANUP:
   if (rval) {
      if (cd) {
         free_colordata(cd);
         free(cd);
      }
      parent_cd->diff_children = (colordata*) NULL;
   }
   if(v1 != v2)
      COLORadjgraph_free(&G);

   return rval;
}

//Initialisiert alle Felder von dsat
int DSATURREC_init(DSATURREC *dsat, int lowerbound)
{
    int rval = 0;
    if (dsat != NULL && lowerbound > 0)
    {
        dsat->color_count = 0;
        dsat->global_lower_bound = lowerbound;
        dsat->uncolored_nodes_number = 0;
	dsat->optimal=0;
	dsat->id_in_child_cd = (int*)NULL;
        dsat->DSAT_value = (int*)NULL;
    }
    else
    {
        rval = 1;
    }
    return rval;
}

//Gibt Knoten und Kanten des COLORproblem sowie die Farbanzahl aus
void print_COLORproblem_edges(colordata *p){
  int i;
  printf("COLORproblem: n: %d, m: %d, Farbanzahl %d\n", p->ncount, p->ecount, p->ccount);
  for(i=0;i<p->ecount;i++){
    printf("Kante %d: %d %d\n",i+1, p->elist[2*i], p->elist[2*i+1]);
  }
  fflush(stdout);
}

//Gibt Knoten und Adjazenzlisten des COLORadjgraph G aus
void print_COLORadjgraph(COLORadjgraph *G)
{
    int i, j;
    printf("COLORadjgraph: \nKnotenanzahl: %d, Kantenanzahl: %d \n", G->ncount, G->ecount);
    for (i = 0; i < G->ncount; i++)
    {
        printf("Die %d Nachbarn des Knoten %d(Farbe: %d): [", G->nodelist[i].degree, i, G->nodelist[i].color);
        for (j = 0; j < G->nodelist[i].degree; j++)
        {
            printf("%d", G->nodelist[i].adj[j]);
            if (j < G->nodelist[i].degree - 1)
                printf(", ");
        }
        printf("]\n");
    }
    fflush(stdout);
}

//Entfernt dominierte Knoten, also v in V, s.d u ex. mit N(v) c= N(u)
//Entfernt Knotenmenge W (aus G) mit |N(w)| < LB und speichert Ergebnis in H = G[V\W]
//node_mapping mapped die Knoten von H und G; i: Knoten in H , node_mapping[i]: Knoten in G
//verwendet, dass Nachbarschaftsliste nache ID aufsteigend sortiert ist
int DSATUR_simplify(COLORadjgraph *G, COLORadjgraph *H, int **node_mapping, int global_lower_bound)
{
    int i, k, j, rval = 0;
    int new_ncount, new_ecount, *new_elist = NULL;
    int *tmp_node_mapping = NULL;
    int dominate;

    //Eingabeparameter überprüfen
    COLOR_IFFREE(*node_mapping, int);

    //alle Knoten(Index im Array) in tmp_node_mapping (mit Wert -1) werden gelöscht, Rest auf Knoten in H gemappt
    tmp_node_mapping = COLOR_SAFE_MALLOC(G->ncount, int);
    COLORcheck_NULL(tmp_node_mapping, "Failed to allocate tmp_node_mapping");

    new_ncount = G->ncount;
    new_ecount = 0;
    //Bestimme alle Knoten mit |N(v)| < LB
    for (i = 0; i < G->ncount; i++)
    {
        if (G->nodelist[i].degree < global_lower_bound)
        {
            tmp_node_mapping[i] = -1;
            new_ncount--;
            //Anzahl der zu löschenden Kanten
            new_ecount += G->nodelist[i].degree;
        }
        else
        {
            tmp_node_mapping[i] = 0;
        }
    }

    //Bestimme alle dominierten Knoten, die nicht eh schon |N(v)|< LB haben, also tmp_node_mapping != -1
    for (i = 0; i < G->ncount; i++)
    {
        if (tmp_node_mapping[i] == 0)
        {
            dominate = 0;
            int neighbor = G->nodelist[i].adj[0]; //erster Nachbar von i(es existiert einer, sonst wäre der Knoten vorher "gelöscht" worden)
            int node;
            //Teste für jeden Nachbarknoten node des neighbor mit größerer Adjazenzliste, ob die von i enthalten ist bis so ein Knoten gefunden wurde
            for (j = 0; j < G->nodelist[neighbor].degree && !dominate; j++)
            {
                node = G->nodelist[neighbor].adj[j];
                if (i != node && tmp_node_mapping[node] != -1 && G->nodelist[i].degree <= G->nodelist[node].degree)
                {
                    //tmp_node_mapping != -1, da es im Fall von -1 einen anderen Knoten mit nicht weniger Kanten gibt, der den Knoten auch dominieren würde

                    //Teste ob Nachbarschaftsliste von i komplett in der von node enthalten ist(verwendet, dass Nachbarschaftsliste sortiert ist)
                    dominate = 1;
                    int l = 0;
                    for (k = 0; (k < G->nodelist[i].degree) && dominate; k++) //geht Nachbarschaftsliste von i durch
                    {                                                         //k ist ein Nachbar von i für den getestet wird, ob er ein Nachbar von node ist
                        while (dominate)                                      //geht Nachbarschaftsliste von node durch
                        {
                            int diff = G->nodelist[i].adj[k] - G->nodelist[node].adj[l];
                            if (diff < 0) //Nachbar ist nicht in der Liste
                            {
                                dominate = 0;
                            }
                            else if (diff > 0) //Nachbar wurde noch nicht in Liste gefunden; Nachbarn sind aufsteigend nach ID sortiert
                            {
                                l++;
                                if (l == G->nodelist[node].degree) //komplette Liste ist durchgelaufen ohne den Nachbarn zu finden
                                    dominate = 0;
                            }
                            else
                            { //Gleichheit, also Knoten haben einen gleichen Nachbarn
                                l++;
                                if (l == G->nodelist[node].degree && k < G->nodelist[i].degree - 1) //komplette Liste von node ist durchgelaufen ohne alle Nachbarn zu finden
                                    dominate = 0;
                                break;
                            }
                        }
                    }
                }
            }

            //Ist es ein dominierter Knoten, dann soll er gelöscht werden
            if (dominate)
            {
                tmp_node_mapping[i] = -1;
                new_ecount += G->nodelist[i].degree; //jede Kante wird doppelt summiert
                new_ncount--;
            }
        }
    }

    //Erstelle mapping mit neuer Knotenliste, sodass später Färbung übertragen werden kann
    (*node_mapping) = COLOR_SAFE_MALLOC(new_ncount, int);
    COLORcheck_NULL((*node_mapping), "Failed to allocate node_mapping");

    //tmp_node_mapping so anpassen, dass Speicherung in H einfacher, node_mapping für Mapping zwischen G und H speichern
    j = 0;
    for (i = 0; i < G->ncount; i++)
    {
        if (tmp_node_mapping[i] != -1)
        {
            tmp_node_mapping[i] = j;  //für leichtere Speicherung in H
            (*node_mapping)[j++] = i; //j ist Knoten in H, der Knoten i in G entspricht

            //andere Enden der zu löschenden Kanten
            for (k = 0; k < G->nodelist[i].degree; k++)
            {
                if (tmp_node_mapping[G->nodelist[i].adj[k]] == -1)
                    new_ecount++; //new_ecount erhöhen, damit wirklich alle Kanten doppelt gezählt wurden
            }
        }
    }

    assert(new_ecount % 2 == 0); //alle zu löschenden Kanten wurden doppelt gezählt, wenn nicht wurde Graph falsch initialisiert
    //richtigen new_ecount ausrechnen als Kanten die bleiben, nicht die die gelöscht werden
    new_ecount = G->ecount - (new_ecount / 2);

    //Neuen Graphen H erstellen
    new_elist = COLOR_SAFE_MALLOC(2 * new_ecount, int);
    COLORcheck_NULL(new_elist, "Failed to allocate new_elist");

    //Speichern der nicht zu löschenden Kanten in new_elist für COLORadjgraph_build
    new_ecount = 0;
    for (i = 0; i < G->ncount; ++i)
    {
        if (tmp_node_mapping[i] != -1) //Knoten, die nicht gelöscht werden sollen
        {
            for (j = 0; j < G->nodelist[i].degree; ++j)
            { //speichere nur Kanten, die nicht zu "zu löschenden Knoten" gehen
                if (G->nodelist[i].adj[j] > i && tmp_node_mapping[G->nodelist[i].adj[j]] != -1)
                { //damit nicht doppelt gespeichert wird, nur Kanten zu "größeren" Adjazenzknoten speichern
                    //Knoten shiften (da einige gelöscht werden) durch tmp_node_mapping(Nachbarn sind immernoch nach ID sortiert)
                    (new_elist)[(new_ecount)*2] = tmp_node_mapping[i];
                    (new_elist)[(new_ecount)*2 + 1] = tmp_node_mapping[G->nodelist[i].adj[j]];
                    (new_ecount)++;
                }
            }
        }
    }

    rval = COLORadjgraph_build(H, new_ncount, new_ecount, new_elist);
    COLORcheck_rval(rval, "COLORadjgraph_build failed");

CLEANUP:
    COLOR_IFFREE(new_elist, int);
    COLOR_IFFREE(tmp_node_mapping, int);
    if (rval)
        COLOR_IFFREE(*node_mapping, int);
    return rval;
}

//Kopiert Farben von G nach sets und speichert Anzahl von color_count nach count
int copy_colors(COLORadjgraph *G, int color_count, COLORset **sets, int *count)
{
    int i, rval = 0;

    COLORcheck_NULL(G, "G is NULL in copy_colors");
    COLORcheck_NULL(sets, "sets is NULL in copy_colors");
    COLORcheck_NULL(count, "count is NULL in copy_colors");
    
    if (COLORdbg_lvl() > 1){
      printf("In copy_colors: Farbanzahl: %d und Farben:\n", color_count);
      for (i = 0; i < G->ncount; i++)
         {
             printf("Knoten %d hat Farbe %d\n", i, G->nodelist[i].color);
         }
      fflush(stdout);
    }
    
    //wenn schon mal Farben gespeichert wurden, wurde schonmal Speicher allokiert und muss freigegeben werden
    if(*count >0){
      COLORfree_sets(sets, count);
    }
    
    if (color_count == 0)
    {
        (*sets) = NULL;
        COLORcheck_NULL(*sets, "color_count is 0 in copy_colors");
    }
    (*sets) = COLOR_SAFE_MALLOC(color_count, COLORset);
    COLORcheck_NULL(*sets, "Failed to allocate sets in copy_colors");

    for (i = 0; i < color_count; i++)
    {
        COLORinit_set(&(*sets)[i]);
    }
    
    //Bestimme Anzahl Knoten in Farbklassen für Speicher allokieren
    for (i = 0; i < G->ncount; i++)
    {
      int c = G->nodelist[i].color;
      if(c > -1){
        ((*sets)[c]).count++;
      }
    }
    
    //Speicher allokieren für Farbklassen
    for (i = 0; i < color_count; i++)
    {
        ((*sets)[i]).members = COLOR_SAFE_MALLOC(((*sets)[i]).count, int);
        COLORcheck_NULL(((*sets)[i]).members, "Failed to allocate sets[i] in copy_colors");

        ((*sets)[i]).count = 0; //wieder auf 0 setzen für Indexierung beim Speichern der members
    }
    
    //alles ist gefärbt, also aktualisiere colorclasses(member + count)
    for (i = 0; i < G->ncount; i++)
    {
        int c = G->nodelist[i].color;
	if(c > -1){
	  ((*sets)[c]).members[((*sets)[c]).count] = i;
	  ((*sets)[c]).count++;
	}
    }
    
    *count = color_count;  
CLEANUP:
    return rval;
}

//Hilfsfunktion für die dsat-Werte der Nachbarn des neu gefärbten Knoten in DSATUR_recursion
//Aktualisiert dsat-Werte der Nachbarn, des neu gefärbten Knoten colored_node und speichert welche Nachbarn aktualisiert wurden in updated_neighbors
//In Stelle colored_node*G->ncount steht Anzahl aktualisierter Nachbarn und danach Indizes dieser Nachbarn in der Nachbarschaftsliste von colored_node
//feasible_colors hat 0, wenn farbe zulässig ist und 1, wenn nicht  
void DSATUR_update_dsat_value_for_neighbors(DSATURREC *dsat, int colored_node, COLORadjgraph *G,int *feasible_colors, int *updated_neighbors)
{
    int i;
    COLORadjnode *node =  &(G->nodelist[colored_node]);
    updated_neighbors[colored_node*G->ncount] = 0;
    for (i = 0; i < node->degree; i++)
    {
        if(!feasible_colors[node->adj[i]*G->ncount+node->color])
        {
            feasible_colors[node->adj[i]*G->ncount+node->color]=1;
	    dsat->DSAT_value[node->adj[i]]++;
	    updated_neighbors[colored_node*G->ncount]++;
            updated_neighbors[colored_node*G->ncount+updated_neighbors[colored_node*G->ncount]]=i;
        }
    }
}

//Hilfsfunktion für DSATUR_recursion um für compute_lower_bound child_cd zu berechnen
//erstellt Graphen mit v1,v2 in derselben oder Farben(v1 und v2 nicht benachbart)
//members in cclasses müssen aufsteigend sortiert sein(gegeben, durch vorherige Verwendung von compute_lower_bound)
int DSATUR_create_branch(colordata *cd, COLORproblem *problem, COLORadjgraph *G, int v1, int v2)
{
    int rval = 0;
    if (!cd->ccount)
        compute_lower_bound(cd, problem);
    
    if(COLORdbg_lvl()> 1){
	printf("In create_branch v1= %d (G: %d)Farbe %d, v2= %d (G: %d) Farbe %d\n",
	       v1,cd->orig_node_ids[v1],G->nodelist[cd->orig_node_ids[v1]].color,v2,cd->orig_node_ids[v2],G->nodelist[cd->orig_node_ids[v2]].color);
	fflush(stdout);
    }
    
    //create differ(modifiziert), wenn v1,v2 unterschiedlichen Farben haben oder gleich sind(gleich um neue Farbe zu erstellen)
    if(G->nodelist[cd->orig_node_ids[v1]].color != G->nodelist[cd->orig_node_ids[v2]].color || v1 == v2)
    {
        rval = DSATUR_create_differ(cd, v1, v2);
        COLORcheck_rval(rval, "Failed in DSATUR_create_differ");

        rval = set_id_and_name(cd->diff_children, problem->ncolordata++, cd->pname);
        COLORcheck_rval(rval, "Failed in set_id_and_name");
    }//create same wenn v1 und v2 dieselbe Farbe haben
    else if (G->nodelist[cd->orig_node_ids[v1]].color == G->nodelist[cd->orig_node_ids[v2]].color)
    {
        rval = create_same(cd, v1, v2);
        COLORcheck_rval(rval, "Failed in create_same");
	
        rval = set_id_and_name(cd->same_children, problem->ncolordata++, cd->pname);
        COLORcheck_rval(rval, "Failure in set_id_and_name");
    }
    
CLEANUP:
    if(rval){
	if(cd->same_children){
	  cd->nsame=1;
	}
	if(cd->diff_children){
	  cd->ndiff=1;
	}
	free_children_data(cd);
    }
    return rval;
}

//Hilfsfunktion zum Sortieren der Knoten in DSATUR_sort_vertices_by_degree
//Sortiert übergebene Knoten (vom Typ sort struct) nach Grad absteigend und Cliqueknoten an erste Stellen
static int DSATUR_comp_by_degree(const void *v1, const void *v2)
{
    const sort node1 = *((const sort *)v1);
    const sort node2 = *((const sort *)v2);
    if ((node1.in_clique && node2.in_clique) || (node1.degree == node2.degree && node1.in_clique == node2.in_clique)) //muss nicht sortiert werden, da beide gleich
        return 0;
    else if (node1.in_clique)
    {
        if (node2.in_clique && node1.degree < node2.degree)
            return 1;
        else
            return -1;
    }
    else
    {
        if (!node2.in_clique && node1.degree > node2.degree)
            return -1;
        else
            return 1;
    }
}

//erstellt einen neuen COLORadjgraph der nach Grad absteigend sortierte Knoten hat und zu vorderst die Cliqueknoten
//node_mapping merkt sich welche Knoten in H (als Indexstelle) dem entsprechenden Knoten in G zugeteilt werden
int DSATUR_sort_vertices_by_degree(COLORadjgraph *G, int **node_mapping, int n, int *clique, int clique_number)
{
    int i, j, k, rval = 0;
    sort *tmp = NULL;
    int *tmp_backwards = NULL;
    int *new_elist = NULL;

    //tmp wird zum Sortieren der Knoten verwendet
    tmp = COLOR_SAFE_MALLOC(G->ncount, sort);
    COLORcheck_NULL(tmp, "Failed to allocate tmp");

    //tmp initialisieren mit bisheriger Sortierung von node_mapping
    for (i = 0; i < G->ncount; i++)
    {
        tmp[i].id = (*node_mapping)[i];
        tmp[i].in_clique = 0;
        tmp[i].degree = G->nodelist[i].degree;
    }

    for (i = 0; i < clique_number; i++)
    {
        tmp[clique[i]].in_clique = 1;
    }

    //Sortiere Knoten mithilfe von DSATUR_comp_by_degree und an erste Stelle stehen die Cliqueknoten
    qsort(tmp, G->ncount, sizeof(sort), DSATUR_comp_by_degree);

    //tmp_backwards wird zum Speichern der Kanten in new_elist gebraucht, um Knoten zu mappen
    tmp_backwards = COLOR_SAFE_MALLOC(n, int);
    COLORcheck_NULL(tmp_backwards, "Failed to allocate tmp_backwards");

    for (i = 0; i < G->ncount; i++)
    {
        tmp_backwards[tmp[i].id] = i;
    }

    //Speichern der Kanten mit den neuen Indizes
    new_elist = COLOR_SAFE_MALLOC(2 * G->ecount, int);
    COLORcheck_NULL(new_elist, "Failed to allocate new_elist");

    k = 0;
    for (i = 0; i < G->ncount; i++)
    {
        for (j = 0; j < G->nodelist[i].degree; j++)
        {
            if (i < G->nodelist[i].adj[j])
            { //existieren keine Schleifen, da Graph vorher schonmal vereinfacht wurde
                new_elist[2 * k] = tmp_backwards[(*node_mapping)[i]];
                new_elist[2 * k + 1] = tmp_backwards[(*node_mapping)[G->nodelist[i].adj[j]]];
                k++;
            }
        }
    }

    //node_mapping auf Sortierung anpassen
    for (i = 0; i < G->ncount; i++)
    {
        (*node_mapping)[i] = tmp[i].id;
    }

    n = G->ncount;
    k = G->ecount;
    COLORadjgraph_free(G);
    rval = COLORadjgraph_build(G, n, k, new_elist);
    COLORcheck_rval(rval, "Failed to build graph");
	
    //Sortiere Nachbarschaftsliste nach ID
    for (i = 0; i < G->ncount; i++)
    {
        qsort(G->nodelist[i].adj, G->nodelist[i].degree, sizeof(int), comp_node_ids);
    }

CLEANUP:
    COLOR_IFFREE(tmp, sort);
    COLOR_IFFREE(tmp_backwards, int);
    COLOR_IFFREE(new_elist, int);
    return rval;
}

//Rekursive Funktion, die von DSATUR_MaxClique aufgerufen wird um maximale Clique zu bestimmen, welche in clique gespeichert wird und max als ihre Größe
//found zeigt an, ob eine Clique gefunden wurde
//clique_sizes_for_nodes[i] speichert Größe der größtmöglichen Clique in {v_i,...,v_n}
//set Knoten müssen nach ID aufsteigend sortiert sein
int DSATUR_MaxClique_recursion(int **clique, int *max, int *found, int *clique_sizes_for_nodes, int *marked_nodes, COLORadjgraph *G, int **set, int *set_size, int size)
{
    int i, j, k, rval = 0;
    int index = 0;
    int cut;
    int node;
    int *newset = (int *)NULL;
    int newset_count;

    //Fehlerbehandlung von Übergabeparametern
    COLORcheck_NULL(clique, "clique is NULL in DSATUR_MaxClique_recursion");
    COLORcheck_NULL(max, "max is NULL in DSATUR_MaxClique_recursion");
    COLORcheck_NULL(found, "found is NULL in DSATUR_MaxClique_recursion");
    COLORcheck_NULL(clique_sizes_for_nodes, "clique_sizes_for_nodes is NULL in DSATUR_MaxClique_recursion");
    COLORcheck_NULL(marked_nodes, "marked_nodes is NULL in DSATUR_MaxClique_recursion");
    COLORcheck_NULL(set_size, "set_size is NULL in DSATUR_MaxClique_recursion");
    //set darf NULL sein, wenn die leere Menge übergeben wurde

    if (*set_size == 0)
    {
        if (size > *max)
        {
            *max = size;

            //Allokiere Speicher neu für gefundene Clique
            COLOR_IFFREE(*clique, int);
            *clique = COLOR_SAFE_MALLOC(*max, int);
            COLORcheck_NULL(*clique, "Failed to allocate clique in DSATUR_MaxClique_recursion");

            //Speichere markierte Knoten in der Clique
            for (i = 0; i < G->ncount; i++)
            {
                if (marked_nodes[i] == 1)
                    (*clique)[index++] = i;
            }
            (*found) = 1;
        }
    }
    else
    {
        //Solange set nicht-leere Menge, also set_size oft, da in jeder Iteration ein Knoten entfernt wird
        for (i = *set_size; i > 0 && size + *set_size > *max; i--)
        {
            //kleinster Knoten im set ist (*set)[0] wegen Sortierung
            node = (*set)[0];
            marked_nodes[node] = 1;

            if (size + clique_sizes_for_nodes[node] <= *max)
            {
                marked_nodes[node] = 0;
                break;
            }

            //der markierte Knoten kann aus set entfernt werden, dafür setze alle Einträge eine Position nach vorne
            for (j = 0; j < *set_size - 1; j++)
            {
                (*set)[j] = (*set)[j + 1];
            }
            (*set_size)--;
            //Lösche letzte Stelle durch reallocaten des Speichers
            if (*set_size)
            {
                (*set) = (int *)realloc((*set), (*set_size) * sizeof(int));
                COLORcheck_NULL((*set), "Failed to reallocate set in DSATUR_MaxClique_recursion");
            }
            else
            {
                COLOR_FREE(*set, int); //Pointer wird auf NULL gesetzt, wenn kein Knoten mehr übrig ist
            }

            // Markieren der Knoten im set, die im Schnitt mit den Nachbarn von node übrig bleiben
            newset_count = *set_size;
            for (j = 0; j < *set_size; j++)
            {
                cut = 0;
                for (k = 0; k < G->nodelist[node].degree; k++)
                {
                    if (G->nodelist[node].adj[k] == (*set)[j])
                    {
                        cut = 1;
                        marked_nodes[(*set)[j]] += 2; //Markierung setzen für Speicherung
                    }
                }
                if (!cut)
                    newset_count--;
            }

            if (newset_count)
            {
                //Entferne node und alle Knoten die nicht Nachbarn von node sind aus set für nächsten Rekursionsschritt
                newset = COLOR_SAFE_MALLOC(newset_count, int);
                COLORcheck_NULL(newset, "Failed to allocate newset in DSATUR_MaxClique_recursion");

                index = 0;
                for (j = 0; j < *set_size; j++)
                {
                    if (marked_nodes[(*set)[j]] > 1)
                    {
                        newset[index++] = (*set)[j];
                        marked_nodes[(*set)[j]] -= 2; //hebe Markierung wieder auf
                    }
                }
            }

            //Rekursiver Aufruf von DSATUR_MaxClique_recursion mit newset
            rval = DSATUR_MaxClique_recursion(clique, max, found, clique_sizes_for_nodes, marked_nodes, G, &newset, &newset_count, size + 1);
            COLORcheck_rval(rval, "Failed to compute DSATUR_MaxClique_recursion");

            COLOR_IFFREE(newset, int);
            marked_nodes[node] = 0;

            if ((*found) == 1)
                break;
        }
    }

CLEANUP:
    if (rval)
        COLOR_IFFREE(newset, int);
    return rval;
}

//Funktion wird nicht mehr verwendet!
//Hilfsfunktion, welche die Größe einer maximalen Clique in lower_bound speichert und Cliqueknoten in clique
//basiert auf Östergards Algorithmus
//sortiert Nachbarschaftslisten nach Id, für schnelle Findung der set-Menge
int DSATUR_MaxClique(int **clique, COLORadjgraph *G, int *lower_bound)
{
    int i, j, found, rval = 0;
    int maximum = 0;
    int set_size;
    int set_index;
    int *clique_sizes_for_nodes = NULL;
    int *set = NULL;
    int *marked_nodes = NULL;
    double starttime = COLORcpu_time();

    COLORcheck_NULL(clique, "clique is NULL in DSATUR_MaxClique");
    COLORcheck_NULL(G, "G is NULL in DSATUR_MaxClique");
    COLORcheck_NULL(lower_bound, "lower_bound is NULL in DSATUR_MaxClique");

    clique_sizes_for_nodes = COLOR_SAFE_MALLOC(G->ncount, int);
    COLORcheck_NULL(clique_sizes_for_nodes, "Failed to allocate clique_size_for_nodes for DSATUR_MaxClique");

    marked_nodes = COLOR_SAFE_MALLOC(G->ncount, int);
    COLORcheck_NULL(marked_nodes, "Failed to allocate marked_nodes for DSATUR_MaxClique");

    //Array für markierte Knoten mit 0 initialisieren
    for (i = 0; i < G->ncount; i++)
    {
        marked_nodes[i] = 0;
    }

    //Sortiere Nachbarschaftsliste nach ID
    for (i = 0; i < G->ncount; i++)
    {
        qsort(G->nodelist[i].adj, G->nodelist[i].degree, sizeof(int), comp_node_ids);
    }

    for (i = G->ncount - 1; i >= 0; i--)
    {
        found = 0;
        marked_nodes[i] = 1;

        //Finden der Schnittmenge zwischen S_i(alle Knoten >=i) und Nachbarn von i 
	//Nachbarn sind aufsteigend nach ID sortiert, daher braucht man nur Index in Nachbarschaftsliste speichern
        set_size = 0;
        for (j = 0; j < G->nodelist[i].degree; j++)
        {
            if (G->nodelist[i].adj[j] >= i)
            {
                set_size = G->nodelist[i].degree - j;
                set_index = j;
                break;
            }
        }

        if (set_size)
        {
            set = COLOR_SAFE_MALLOC(set_size, int);
            COLORcheck_NULL(set, "Failed to allocate set for DSATUR_MaxClique");

            for (j = set_index; j < G->nodelist[i].degree; j++)
            {
                set[j - (set_index)] = G->nodelist[i].adj[j];
            }
        }
        else
        {
            set = (int *)NULL;
        }

        //Aufruf von DSATUR_MaxClique_recursion um Clique zu bestimmen
        rval = DSATUR_MaxClique_recursion(clique, &maximum, &found, clique_sizes_for_nodes, marked_nodes, G, &set, &set_size, 1);
        COLORcheck_rval(rval, "Failed to compute DSATUR_MaxClique_recursion");

        //Speicher von set wieder freigeben für nächsten Schleifendurchlauf
        COLOR_IFFREE(set, int);

        marked_nodes[i] = 0; //Für nächsten Schleifendurchlauf Markierung entfernen

        clique_sizes_for_nodes[i] = maximum;
    }

    //Lower bound aktualisieren
    *lower_bound = *lower_bound > maximum ? *lower_bound : maximum;

CLEANUP:

    COLOR_IFFREE(marked_nodes, int);
    COLOR_IFFREE(clique_sizes_for_nodes, int);
    if (rval)
        COLOR_IFFREE(set, int);
    printf("Zeit von Maxclique: %lf\n", COLORcpu_time() - starttime);
    return rval;
}

//wählt und berechnet geeignete globale untere Schranke (LB) vor Beginn der Rekursion
//Nimmt für Graphen der Dichte (nach Vereinfachung) <= 0.8 Heuristik basierend auf Östergards Algorithmus und sonst fraktionale chromatische Zahl
//Speichert Cliqueknoten (bei compute_lower_bound 2-er Clique) in clique
int DSATUR_choose_and_compute_lower_bound(DSATURREC *dsat,int **clique, COLORadjgraph *G, COLORproblem *p)
{
    int rval=0;
    int* elistc=NULL;
    COLORset* newsets  = (COLORset*) NULL;
    COLORNWT *nweights = (int *) NULL;
    int nnewsets =0;
    double starttime = COLORcpu_time();
    double density = ((double) G->ecount) * 2.0 / ((double) (G->ncount * (G->ncount - 1)));

/*	Wird nicht mehr verwendet, da Heuristik im Allgemeinen schneller als DSATUR_MaxClique ist
 *    if(density < 0.0){
      //MaxCLiquealgorithmus nach Östergard für COLORadjgraph
      rval = DSATUR_MaxClique(clique, G, &(dsat->global_lower_bound));
      COLORcheck_rval(rval, "Failed to DSATUR_MaxClique");
      //color_count und uncolored_nodes_number initialisieren
      dsat->color_count = dsat->global_lower_bound;
      printf("Zeit: %lf für lower_bound: %d mit MaxClique\n", COLORcpu_time()-starttime, dsat->global_lower_bound);
    }*/
    if(density < 0.8){
        //Benutze COLORclique_ostergard als Heuristik
        int ecountc;
	int i;
	int solutionweight = 1;
	
	nweights = COLOR_SAFE_MALLOC (G->ncount, int);
        COLORcheck_NULL (nweights, "out of memory for nweights");
        for (i = 0; i < G->ncount; ++i) {
          nweights[i] = 1;
        }

        COLORadjgraph_extract_edgelist(&ecountc,&elistc,G);
	
	rval = COLORclique_ostergard(&newsets, &nnewsets, G->ncount, G->ecount, elistc, nweights, COLOR_MAXINT, &solutionweight, 100);
	COLORcheck_rval(rval,"COLORclique_ostergard failed");
	
	if(newsets->count < 2){
	  *clique = COLOR_SAFE_MALLOC(2, int);
	  COLORcheck_NULL(*clique, "no memory allocated for clique");

	  (*clique)[0] = 0;
	  (*clique)[1] = G->nodelist[0].adj[0];

	  dsat->color_count = 2;
	}
	else{
	  dsat->color_count = newsets->count;
	  *clique = COLOR_SAFE_MALLOC (newsets->count, int);
	  COLORcheck_NULL (*clique, "out of memory for clique");
	  
	  for(i=0;i<dsat->color_count;i++)
	  {
	    (*clique)[i]=newsets->members[i];
	  }
	}
	dsat->global_lower_bound = p->root_cd.lower_bound > dsat->color_count ? p->root_cd.lower_bound : dsat->color_count;
	p->root_cd.lower_bound = dsat->global_lower_bound;
	printf("Time for lower_bound: %lf seconds with COLORclique_ostergard and lower_bound %d\n", COLORcpu_time()-starttime, dsat->color_count);
    } else{
        //Benutze compute_lower_bound
	//compute_lower_bound berechnet die LB auf dem noch nicht vereinfachten Graphen
        rval = compute_lower_bound(&(p->root_cd),p);
        COLORcheck_rval(rval,"Failed to compute_lower_bound");
	
        dsat->global_lower_bound = p->root_cd.lower_bound;

        *clique = COLOR_SAFE_MALLOC(2, int);
        COLORcheck_NULL(*clique, "no memory allocated for clique");

        (*clique)[0] = 0;
        (*clique)[1] = G->nodelist[0].adj[0];

        dsat->color_count = 2;
	printf("Time for lower_bound: %lf seconds with compute_lower_bound and lower_bound %d\n", COLORcpu_time()-starttime, dsat->color_count);
    }
    dsat->uncolored_nodes_number = G->ncount - dsat->color_count;
CLEANUP:
    COLOR_IFFREE(nweights, COLORNWT);
    COLORfree_sets(&newsets, &nnewsets);
    COLOR_IFFREE(elistc, int);
    return rval;
}

//Erstellt neues COLORproblem für die Rekursion (übernimmt einige Felder von problem)
//Für jeden Cliqueknoten wird ein stable set gespeichert und eine erste Färbung initialisiert
int DSATUR_new_COLORproblem(COLORproblem *new_problem, COLORproblem *problem, COLORadjgraph *G, int color_count)
{
    int i, rval = 0;
    colordata *cd = &(new_problem->root_cd);
    colordata *old_cd = &(problem->root_cd);

    cd->id = 0;
    new_problem->ncolordata = 1;
    new_problem->global_upper_bound = problem->global_upper_bound;

    cd->ncount = G->ncount;
    rval = COLORadjgraph_extract_edgelist(&cd->ecount, &cd->elist, G);
    COLORcheck_rval(rval, "Failed in COLORadjgraph_extract_edgelist");

    cd->orig_node_ids = (int *)COLOR_SAFE_MALLOC(cd->ncount, int);
    COLORcheck_NULL(cd->orig_node_ids, "Failed to allocate cd->orig_node_ids");
    for (i = 0; i < cd->ncount; i++)
    {
        cd->orig_node_ids[i] = i;
    }

    cd->upper_bound = old_cd->upper_bound;
    cd->lower_bound = old_cd->lower_bound;
    
    //speichern der cliquesets in cclasses
    cd->cclasses = COLOR_SAFE_MALLOC(color_count, COLORset);
    COLORcheck_NULL(cd->cclasses, "Failed to allocate cclasses");
    cd->ccount = color_count;
    cd->gallocated = cd->ccount;

    for (i = 0; i < color_count; i++)
    {
      COLORinit_set(&(cd->cclasses)[i]);
      cd->cclasses[i].count = 1;
      (cd->cclasses[i]).members = COLOR_SAFE_MALLOC(1, int);
      COLORcheck_NULL((cd->cclasses[i]).members, "Failed to allocate sets[i] in DSATUR_new_COLORproblem");
      cd->cclasses[i].members[0] = i;
    }
	
    rval = COLORdsatur (cd->ncount, cd->ecount, cd->elist, &(cd->ccount), &(cd->cclasses));
    COLORcheck_rval (rval, "COLORdsatur failed");
    COLORcopy_sets(&(cd->bestcolors),&(cd->nbestcolors), cd->cclasses,cd->ccount);
    COLORcheck_coloring(cd->bestcolors,cd->nbestcolors, cd->ncount, cd->ecount, cd->elist);
    
    cd->upper_bound = cd->nbestcolors < cd->upper_bound ? cd->nbestcolors : cd->upper_bound;
    new_problem->global_upper_bound = cd->upper_bound;
    cd->gallocated = cd->ccount;

    new_problem->key_mult = (COLORNWT_MAX - 1) / cd->ncount;

CLEANUP:
    return rval;
}

//Algorithmus 1:
//eigene DSATURfunktion mit fraktionaler chromatischer Zahl als Grenze
//speichert optimale Färbung in ncolors und colorclasses
//erhält Graphen durch COLORproblem
//berechnet zuerst vereinfachten Graphen H und berechnet mithilfe von DSATUR_choose_and_compute_lower_bound eine globale untere Schranke
//sortiert dann Knoten und ruft dann erst DSATUR_recursion auf
//am Ende wird Färbung auf G übertragen und in colorclasses gespeichert
int DSATUR(COLORproblem *problem, int *ncolors, COLORset **colorclasses)
{
    int i, j, rval = 0;
    colordata *cd = &(problem->root_cd);
    COLORadjgraph G;
    COLORadjgraph H; //H ist Hilfsgraph, der durch Vereinfachung von G entsteht
    DSATURREC dsat;
    COLORproblem new_problem;
    COLORproblem_init(&new_problem);
    colordata *ncd = &(new_problem.root_cd);
    int *node_mapping = (int *)NULL;
    int **clique = (int **)NULL;
    int *feasible_colors = (int *)NULL;
    int *colored_node_updated_neighbors = NULL;
    double *lbzeit = NULL;
    int *aufrufanzahl = NULL;

    COLORcheck_NULL(problem, "problem is NULL in DSATUR");
    COLORcheck_NULL(ncolors, "ncolors is NULL in DSATUR");
    COLORcheck_NULL(colorclasses, "colorclasses is NULL in DSATUR");

    clique = COLOR_SAFE_MALLOC(1, int *);
    COLORcheck_NULL(clique, "Failed to allocate clique");
    *clique = (int *)NULL;

    //Adjazenzgraph G bauen
    rval = COLORadjgraph_build(&G, cd->ncount, cd->ecount, cd->elist);
    COLORcheck_rval(rval, "Failed to build_graph");

    //Graph vereinfachen: 1. Mehrfachkanten & Schleifen entfernen (sortiert auch Nachbarschaftslisten nach aufsteigender ID)
    rval = COLORadjgraph_simplify(&G);
    COLORcheck_rval(rval, "Failed to simplify graph");

    //Farben initialisieren in G
    for (i = 0; i < G.ncount; i++)
        G.nodelist[i].color = -1;

    if (COLORdbg_lvl() > 0)
    {
        printf("G in DSATUR ");
        print_COLORadjgraph(&G);
        fflush(stdout);
    }

    //Upper Bound initialisieren (fängt auch ab, wenn UB > n ist)
    problem->global_upper_bound = cd->upper_bound < cd->ncount ? cd->upper_bound : cd->ncount;

    rval = DSATURREC_init(&dsat, cd->lower_bound);
    COLORcheck_rval(rval, "Failed to init DSATURREC");

    //G vereinfachen und Ergebnis in H speichern: 2. Knoten v mit |N(v)| < LB und dominierte Knoten löschen
    //braucht, dass Nachbarschaftsliste nach ID aufsteigend sortiert ist(hier gegeben durch COLORadjgraph_simplify)
    rval = DSATUR_simplify(&G, &H, &node_mapping, cd->lower_bound);
    COLORcheck_rval(rval, "Failed in DSATUR_simplify");

    //Wählt geeignete untere Schranke und bestimmt mindestens eine 2-er clique
    rval = DSATUR_choose_and_compute_lower_bound(&dsat, clique, &H, problem);
    COLORcheck_rval(rval, "Failed in DSATUR_choose_and_compute_lower_bound");

    //Arrays in dsat initialisieren, die H.ncount brauchen
    dsat.DSAT_value = COLOR_SAFE_MALLOC(H.ncount, int);
    COLORcheck_NULL(dsat.DSAT_value, "Failed to allocate DSAT_value");
    dsat.id_in_child_cd = COLOR_SAFE_MALLOC(H.ncount, int);
    COLORcheck_NULL(dsat.id_in_child_cd, "Failed to allocate id_in_child_cd");
    
    for(i=0;i<H.ncount;i++){
      dsat.DSAT_value[i]=0;
      dsat.id_in_child_cd[i]=-1;
    }

    //G vereinfachen: 3. Knoten sortieren: erst nach Clique und dann nach Grad absteigend
    rval = DSATUR_sort_vertices_by_degree(&H, &node_mapping, G.ncount, *clique, dsat.color_count);
    COLORcheck_rval(rval, "Failed in DSATUR_sort_vertices_by_degree");

    //wird mit 0 initialisiert, da anfangs jede Farbe für jeden Knoten zulässig ist
    //speichert an stelle i*H.ncount+j 0 wenn Farbe j für Knoten i zulässig
    feasible_colors = (int *)calloc(H.ncount * H.ncount, sizeof(int));
    COLORcheck_NULL(feasible_colors, "Failed to allocate feasible_colors");

    //Cliqueknoten, welche nach der Sortierung die ersten Knoten in H sind, färben und für debug ausgeben und Hilfsarrays aktualisieren
    for (i = 0; i < dsat.color_count; i++)
    {
        H.nodelist[i].color = i;
	dsat.id_in_child_cd[i] = i;
        for (j = 0; j < H.nodelist[i].degree; j++)
        {
            feasible_colors[H.nodelist[i].adj[j] * H.ncount + i] = 1;
            dsat.DSAT_value[H.nodelist[i].adj[j]]++;
        }
    }
    //Farben auf -1 setzen in H
    for (i = dsat.color_count; i < H.ncount; i++)
    {
        H.nodelist[i].color = -1;
    }
    
    if (COLORdbg_lvl() > 0)
    {
        printf("%d Cliqueknoten:\t", dsat.color_count);
        for (i = 0; i < dsat.color_count; i++)
            printf("%d ", (*clique)[i]);
        printf("\n");

        printf("H in DSATUR ");
        print_COLORadjgraph(&H);
        fflush(stdout);
    }

    //colored_node_updated_neighbors speichert welche Nachbarn von colored_node einen aktualisierten dsat-Wert bekommen
    //Speichert an Stelle i*H.ncount die Anzahl der Nachbarn von i, die aktualisiert wurden und danach die entsprechenden Nachbarn
    colored_node_updated_neighbors = COLOR_SAFE_MALLOC(H.ncount * H.ncount, int);
    COLORcheck_NULL(colored_node_updated_neighbors, "Failed to allocate colored_node_updated_neighbors for DSATUR_recursion");

    //erstellt neues COLORproblem, das den Graph H abbildet
    rval = DSATUR_new_COLORproblem(&new_problem, problem, &H, dsat.color_count);
    COLORcheck_rval(rval, "Failed in DSATUR_new_COLORproblem");

    //Variablen zum Auslesen der Zeit und Aufrufanzahl
    lbzeit = (double *)calloc(1, sizeof(double));
    aufrufanzahl = (int *)malloc(sizeof(int));
    *aufrufanzahl = 1;
    double dsaturtime = COLORcpu_time();
    
    assert(dsat.color_count > 1);
    //Aufruf der rekursiven DSATUR-Funktion
    rval = DSATUR_recursion(ncd, &new_problem, &H, &dsat, feasible_colors, aufrufanzahl, lbzeit, colored_node_updated_neighbors);
    COLORcheck_rval(rval, "Failed in DSATUR_recursion");

    printf("Time for DSATUR_recursion: %lf seconds of which time for compute_lower_bound: %lf seconds\n", COLORcpu_time() - dsaturtime, *lbzeit);
    nodes += *aufrufanzahl;
    
    //in H gelöschte Knoten => jetzt färben in G
    //erst bereits Farben vorhandener Knoten kopieren, dazu wird node_mapping verwendet
    for (i = 0; i < ncd->nbestcolors; i++)
    {
        for (j = 0; j < ncd->bestcolors[i].count; j++)
        {
            G.nodelist[node_mapping[ncd->bestcolors[i].members[j]]].color = i;
        }
    }
    //Rest in G färben, wenn überhaupt Knoten gelöscht wurden durch die Vereinfachung
    if (G.ncount > H.ncount)
    {
        for (i = 0; i < G.ncount; i++)
        { //Ändert an chromatischer Zahl nichts, da nur "leicht zu färbende" Knoten gelöscht wurden
            if (G.nodelist[i].color == -1)
                color_node(&G, i);
        }
    }
    
    //Speichere Farben neu in colorclasses und ncolors, damit id der Knoten von G und nicht von H verwendet werden
    rval = copy_colors(&G, ncd->nbestcolors, colorclasses, ncolors);
    COLORcheck_rval(rval, "Failed to copy_colors in DSATUR");

    //Für Ausgabe in main aktualisiere cd->upper_bound und setze cd->lower_bound auf die globale untere Schranke
    cd->upper_bound = new_problem.global_upper_bound;
    cd->lower_bound = dsat.global_lower_bound;

    rval = COLORcheck_coloring(*colorclasses, *ncolors, cd->ncount, cd->ecount, cd->elist);
    COLORcheck_rval(rval, "Failed to verify coloring in DSATUR");
	
    if (COLORdbg_lvl() > 0)
    {
        printf("Farben in DSATUR G:\n");
        for (i = 0; i < G.ncount; i++)
        {
            printf("Knoten %d hat Farbe %d\n", i, G.nodelist[i].color);
        }
        fflush(stdout);
    }

CLEANUP:
    COLORproblem_free(&new_problem);
    COLORfree_sets(&cd->bestcolors, &cd->nbestcolors);
    COLOR_IFFREE(*clique, int);
    COLOR_IFFREE(clique, int *);
    COLORadjgraph_free(&G);
    COLORadjgraph_free(&H);
    COLOR_IFFREE(node_mapping, int);
    COLOR_IFFREE(colored_node_updated_neighbors, int);
    COLOR_IFFREE(feasible_colors, int);
    COLOR_IFFREE(dsat.DSAT_value, int);
    COLOR_IFFREE(dsat.id_in_child_cd, int);
    COLOR_IFFREE(aufrufanzahl, int);
    COLOR_IFFREE(lbzeit, double);

    return rval;
}

//Algorithmus 2:
//berechnet rekursiv optimale Färbung und speichert diese in problem->root_cd.bestcolors und nbestcolors
//speichern von UB, LB, bestcolors in problem und Berechnen der Farben etc. mit G
//erstellt child_cd, dass aktuelle partielle Färbung darstellt für compute_lower_bound mithilfe von id_in_child_cd und create_branch
int DSATUR_recursion(colordata *child_cd, COLORproblem *problem, COLORadjgraph *G, DSATURREC *dsat, int *feasible_colors, int *aufrufanzahl, double *lbzeit, int *colored_node_updated_neighbors)
{
    int i, j, rval = 0;
    colordata *cd = &(problem->root_cd);
    int local_lower_bound;
    int colored_node = 0;
    
    //Sind alle Knoten gefärbt?
    if (dsat->uncolored_nodes_number == 0)
    {
        if (dsat->color_count < problem->global_upper_bound)
        {
            rval = copy_colors(G, dsat->color_count, &cd->bestcolors, &cd->nbestcolors);
            COLORcheck_rval(rval, "Failed to copy colors");
            problem->global_upper_bound = dsat->color_count;
            if (cd->nbestcolors == dsat->global_lower_bound) //optimale Lösung wurde gefunden
            {
                dsat->optimal = 1;
            }
        }
    }
    else //rekursiv einen Knoten färben
    {
        local_lower_bound = dsat->global_lower_bound > dsat->color_count ? dsat->global_lower_bound : dsat->color_count;
        child_cd->lower_bound = child_cd->lower_bound > local_lower_bound? child_cd->lower_bound : local_lower_bound;
	
        if (child_cd->parent)
        { //Berechne neue lokale untere Schranke mit compute_lower_bound

	    //build_lp damit optimize benutzt werden kann
	    if (!child_cd->lp)
	    {
		rval = build_lp(child_cd);
		COLORcheck_rval(rval, "build_lp failed");
	    }

	    rval = COLORlp_optimize(child_cd->lp);
	    COLORcheck_rval(rval, "COLORlp_optimize failed");
	  
            double starttime = COLORcpu_time();
            rval = compute_lower_bound(child_cd, problem);
            COLORcheck_rval(rval, "Failed to compute lower bound");
            *lbzeit += COLORcpu_time() - starttime;

            if (COLORdbg_lvl() > 1)
            {
                printf("Zeit für compute_lower_bound %lf\n", COLORcpu_time() - starttime);
                printf("Für %d verwendete Farben berechnete fraktionale chromatische Zahl: %d\n", dsat->color_count, child_cd->lower_bound);
		fflush(stdout);
            }
            local_lower_bound = local_lower_bound > child_cd->lower_bound ? local_lower_bound : child_cd->lower_bound;

        }
        //Färbung kann vielleicht jetzt schon verworfen werden
        if (local_lower_bound >= problem->global_upper_bound)
        {
            goto CLEANUP;
        }

        //Berechneung des zu färbenden Knoten colored_node
        //Bestimme erstmal einen ungefärbten Knoten, um auch maxsatdegree initialisieren zu können
        colored_node = 0;
        if (G->nodelist[colored_node].color > -1)
        {
            for (i = 1; i < G->ncount; i++)
            {
                if (G->nodelist[i].color == -1)
                {
                    colored_node = i;
                    break;
                }
            }
        }
        int maxsatdegree = dsat->DSAT_value[colored_node];

        //select uncolored vertex mit maximalen DSAT_value
        for (i = colored_node + 1; i < G->ncount; i++)
        { //tiebreaking mit maxdegree; DSAT_value wurde mit 0 vorinitialisiert
            if (G->nodelist[i].color == -1 && (dsat->DSAT_value[i] > maxsatdegree || (dsat->DSAT_value[i] == maxsatdegree && G->nodelist[i].degree > G->nodelist[colored_node].degree)))
            {
                colored_node = i;
                maxsatdegree = dsat->DSAT_value[i];
            }
        }
        
        //Bestimme id von colored_node in child_cd; child_cd hat vielleicht weniger Knoten
        int colored_node_in_cd = colored_node < child_cd->ncount ? colored_node : child_cd->ncount - 1;
        while (child_cd->orig_node_ids[colored_node_in_cd] != colored_node)
        {
            colored_node_in_cd--;
        }
        //Speichert für neue Farbe in was_colored die letzte Farbe, in der colored_node gefärbt war 
        int was_colored = -1;
        int old_color_count = dsat->color_count;
        //Färbe colored_node mit jeder zulässigen Farbe
        for (i = 0; i <= old_color_count; ++i)
        {
            if (!feasible_colors[colored_node * G->ncount + i])
            {
                G->nodelist[colored_node].color = i;
                if (i == dsat->color_count)
                { //wenn es die neue Farbe ist
                    dsat->color_count++;
                }
                //Wenn die obere Schranke nicht überschritten wird, dann gehe tiefer in die Rekursion
                if (dsat->color_count < problem->global_upper_bound)
                {
                    dsat->uncolored_nodes_number--;
		    //Aktualisiere dsat-Werte der Nachbarn von colored_node
                    DSATUR_update_dsat_value_for_neighbors(dsat, colored_node, G, feasible_colors, colored_node_updated_neighbors);
		    
		    *aufrufanzahl += 1;
		    int node_of_old_color = -1;
		    int tmp = -1;
                    //erstelle für compute_lower_bound den farbkomprimierten Graphen
                    if (dsat->uncolored_nodes_number) //sonst wird im nächsten schritt nur die Färbung gespeichert
                    {
                        //Wenn es die neue Farbe ist, kann ein nicht benachbarter gefärbter Knoten gewählt werden
                        if (dsat->color_count > old_color_count)
                        {
			  if(was_colored > -1){//was_colored ist die nicht-benachbarte Farbe
			    //Bestimme Knoten, der entsprechend gefärbt ist, in child_cd
			    node_of_old_color = dsat->id_in_child_cd[was_colored];
			  }
			  else{
			  //Wenn es keine unbenachbarte Farbe gibt, muss die neue Farbe hinzugefügt werden
			  //dafür benutze DSATUR_create_differ mit v1=v2
			    node_of_old_color = colored_node_in_cd;
			  }
			  //id_in_child_cd gibt die id der Farbe in child_cd an
			  dsat->id_in_child_cd[i]= colored_node_in_cd;
                        }
                        else
                        {
			  //wähle Knoten aus Farbklasse der gewählten Farbe
			    node_of_old_color = dsat->id_in_child_cd[i];
			    tmp = node_of_old_color;
			    //v1>v2 erzeugt Probleme beim shiften der Indices in create_same
			    if(node_of_old_color > colored_node_in_cd){//Tausche die Knoten
			      node_of_old_color = colored_node_in_cd;
			      colored_node_in_cd = tmp;
			    }
                        }
                        rval = DSATUR_create_branch(child_cd, problem, G, node_of_old_color, colored_node_in_cd);
                        COLORcheck_rval(rval, "Failed in DSATUR_create_branch");
                    }
                    //Aktualisiere Datenstruktur für nächsten Rekursionsaufruf und setze es danach wie vorher
                    if (child_cd->nsame)
                    {//Fall 1: Farbe existierte schon
			dsat->id_in_child_cd[i] = node_of_old_color;
			for(j=0;j<old_color_count;j++){
			  if(dsat->id_in_child_cd[j] > colored_node_in_cd){
			    dsat->id_in_child_cd[j]--;
			  }
			}
		      
                        rval = DSATUR_recursion(child_cd->same_children, problem, G, dsat, feasible_colors, aufrufanzahl, lbzeit, colored_node_updated_neighbors);
			
			for(j=0;j<old_color_count;j++){
			  if(dsat->id_in_child_cd[j] >= colored_node_in_cd){
			    dsat->id_in_child_cd[j]++;
			  }
			}
			dsat->id_in_child_cd[i] = tmp;
			colored_node_in_cd = node_of_old_color == tmp? colored_node_in_cd : node_of_old_color;
		    }else if(child_cd->ndiff){//Fall 2: neue Farbe wurde benutzt
                        rval = DSATUR_recursion(child_cd->diff_children, problem, G, dsat, feasible_colors, aufrufanzahl, lbzeit, colored_node_updated_neighbors);
                    
			dsat->id_in_child_cd[i] = -1;
		    }
                    else
                    {//Fall 3: DSATUR_create_branch wurde nicht aufgerufen, da im nächsten Schritt nur die Färbung gespeichert wird
                        rval = DSATUR_recursion(child_cd, problem, G, dsat, feasible_colors, aufrufanzahl, lbzeit, colored_node_updated_neighbors);
                    }
                    COLORcheck_rval(rval, "Failed in DSATUR_recursion");
		    
		    //Damit in der nächsten Iteration der for-Schleife wieder ein child generiert werden kann
                    free_children_data(child_cd);
		    
                    if (dsat->optimal == 1)
                        return rval;
		    
                    //Knoten wieder "entfärben"
                    dsat->uncolored_nodes_number++;
                    //dsat_value der Nachbarn wieder zurücksetzen mithilfe von colored_node_updated_neighbors
                    for (j = 1; j < colored_node_updated_neighbors[colored_node * G->ncount] + 1; j++)
                    { //wer geupdated wurde bekam eine Farbe mehr, die er vorher in der Nachbarschaft nicht hatte
                        int neighbor = G->nodelist[colored_node].adj[colored_node_updated_neighbors[colored_node * G->ncount + j]];
                        dsat->DSAT_value[neighbor]--;
                        feasible_colors[neighbor * G->ncount + i] = 0;
                    }
                }
                G->nodelist[colored_node].color = -1;
		was_colored = i;
            }
        }
        dsat->color_count--;
    }

CLEANUP:
    if(rval)free_children_data(child_cd);
    return rval;
}

//main die nur Graphen einliest und dann DSATUR aufruft und dann Speicher freigibt
int main(int ac, char **av)
{
    int rval = 0;
    double start_time = COLORcpu_time();
    double tot_rtime;
    
    COLORproblem colorproblem;
    COLORparms *parms = &(colorproblem.parms);
    colordata *cd = &(colorproblem.root_cd);

    int ncolors = 0;
    COLORset *colorclasses = (COLORset *)NULL;

    rval = COLORprogram_header(ac, av);
    COLORcheck_rval(rval, "Failed in COLORprogram_header");

    COLORproblem_init(&colorproblem);
    cd->id = 0;
    colorproblem.ncolordata = 1;

    rval = parseargs(ac, av, parms);
    if (rval)
        goto CLEANUP;

    //setzt durch Kommandozeilenparameter übergebene upper_bound oder Default
    cd->upper_bound = parms->initial_upper_bound;
    get_problem_name(cd->pname, parms->edgefile);

    printf("Rounding strategy: ");
    switch (parms->rounding_strategy)
    {
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

    if (COLORdbg_lvl() > 1)
        printf("Debugging turned on\n");
    fflush(stdout);

    rval = COLORread_dimacs(parms->edgefile, &(cd->ncount), &(cd->ecount),
                            &(cd->elist), (int **)NULL);
    COLORcheck_rval(rval, "COLORread_dimacs failed");

    if (cd->upper_bound > cd->ncount)
    {
        cd->upper_bound = cd->ncount;
    }

    if (colorproblem.parms.backupdir)
    {
        recover_colordata(cd, &colorproblem);
    }
    
    rval = DSATUR(&colorproblem, &ncolors, &colorclasses);
    COLORcheck_rval(rval, "Failed in DSATUR");
    
    if (cd->nbestcolors == cd->upper_bound)
    {
        printf("Opt Colors: %d\n", cd->nbestcolors);
        fflush(stdout);
        print_colors(colorclasses, ncolors);
    }
    else if (cd->lower_bound == parms->initial_upper_bound)
    {
        printf("Lower bound reached predefined upper bound %d.\n",
               parms->initial_upper_bound);
    }
    else
    {
        printf("Finished with LB %d and UB %d.\n",
               cd->lower_bound, cd->upper_bound);
    }
    tot_rtime = COLORcpu_time() - start_time;
    printf("Computing coloring took %f seconds. colorcount = %d and nodes = %ld\n", tot_rtime, ncolors, nodes);
    
CLEANUP:

    COLORproblem_free(&colorproblem);
    COLORfree_sets(&colorclasses, &ncolors);

    return rval;
}
