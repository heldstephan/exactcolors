#ifndef __GRAPH_H
#define __GRAPH_H


typedef struct node {
    int *adj;
    int degree;
    int color;
} node;

typedef struct graph {
    node *nodelist;
    int  *adjspace;
    int ncount;
    int ecount;
} graph;


int  COLORadjgraph_build(graph* G,int ncount,int ecount, const int elist[]);
int  COLORadjgraph_copy(graph* Gdst, const graph* Gsrc);
int  COLORadjgraph_delete_unweighted(graph* G, int** new_nweights,const int nweights[]);
int  COLORadjgraph_build_complement(graph* Gc, const graph* G);
void COLORadjgraph_free(graph* G);
int  COLORadjgraph_simplify(graph* G);
int  COLORadjgraph_extract_edgelist(int* ecount, int* elist[], const graph* G);
void COLORadjgraph_sort_adjlists_by_id(graph* G);


#endif
