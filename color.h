#ifndef __COLOR_H
#define __COLOR_H
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

#include "color_defs.h"

typedef struct COLORset {
    int count;
    int age;
    int *members;
    struct COLORset *next;
} COLORset;

#define COLOR_PRANDMAX 1000000007

typedef struct COLORrandstate {
    int a;
    int b;
    int arr[55];
} COLORrandstate;


int COLORdbg_lvl(void);

int COLORgreedy (int ncount, int ecount, int *elist, int *ncolors,
        COLORset **colorclasses);
int COLORclique_enum (COLORset** newsets, int *nnewsets, int ncount,
        int ecount, int *elist, int *weights, int cutoff, int *pval);
int COLORclique_ostergard (COLORset **newsets, int *nnewsets, int ncount,
        int ecount, int *elist, int *weights, int cutoff, int *pval);

void COLORinit_set (COLORset *s);
void COLORfree_set (COLORset *s);
void COLORfree_sets (COLORset **s,int* nsets);
int  COLORcopy_sets (COLORset **dsts,int* nsets,
                     const COLORset *src_s, int src_nsets);
void COLORunique_sets (COLORset **s,int* nsets);
int  COLORcheck_set(COLORset* set, int ncount, int ecount, const int elist[]);

/** Test whether (set,ncount) defines a feasible coloring for (ncount,elist,ecount).*/ 
int  COLORcheck_coloring(COLORset* set, int ccount, int ncount, int ecount, const int elist[]);

void COLORutil_sprand (int seed, COLORrandstate *r);
int COLORutil_lprand (COLORrandstate *r);
double COLORutil_zeit (void);
void COLORutil_quicksort (int *len, int n);
void COLORutil_quicksort_reverse (int *len, int n);
void COLORutil_perm_quicksort (int *perm, int *len, int n);

void COLORset_quicksort (COLORset *cclasses, int ccount);



#endif  /* __COLOR_H */
