#ifndef __COLOR_H
#define __COLOR_H

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

void *COLORutil_allocrus (size_t size);
void COLORutil_freerus (void *p);
void COLORinit_set (COLORset *s);
void COLORfree_set (COLORset *s);
void COLORfree_sets (COLORset **s,int* nsets);
int  COLORcheck_set(COLORset* set, int ncount, int ecount, const int elist[]);
void COLORutil_sprand (int seed, COLORrandstate *r);
int COLORutil_lprand (COLORrandstate *r);
double COLORutil_zeit (void);
void COLORutil_quicksort (int *len, int n);
void COLORutil_quicksort_reverse (int *len, int n);
void COLORutil_perm_quicksort (int *perm, int *len, int n);



#endif  /* __COLOR_H */
