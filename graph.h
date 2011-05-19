#ifndef __GRAPH_H
#define __GRAPH_H
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

typedef struct COLORadjnode {
    int* adj;
    int  degree;
    int  color;
} COLORadjnode;

typedef struct COLORadjgraph {
    COLORadjnode* nodelist;
    int*          adjspace;
    int           ncount;
    int           ecount;
} COLORadjgraph;


int  COLORadjgraph_build(COLORadjgraph* G,int ncount,int ecount, const int elist[]);
int  COLORadjgraph_copy(COLORadjgraph* Gdst, const COLORadjgraph* Gsrc);
int  COLORadjgraph_delete_unweighted(COLORadjgraph* G, int** new_nweights,
                                     const int nweights[]);
int  COLORadjgraph_build_complement(COLORadjgraph* Gc, const COLORadjgraph* G);
void COLORadjgraph_init(COLORadjgraph* G);
void COLORadjgraph_free(COLORadjgraph* G);
int  COLORadjgraph_simplify(COLORadjgraph* G);
int  COLORadjgraph_extract_edgelist(int* ecount, int* elist[], const COLORadjgraph* G);
void COLORadjgraph_sort_adjlists_by_id(COLORadjgraph* G);
int  COLORread_dimacs (char *f, int *pncount, int *pecount, int **pelist,
                       int **pnweights);
int  COLORedge_stat(const COLORadjgraph* G);

int  COLORgraph_print(int ecount, const int elist[]);

int  COLORcheck_connectedness(const COLORadjgraph* G);
 
#endif
