#ifndef __COLOR_H
#define __COLOR_H

typedef struct COLORset {
    int count;
    int *members;
} COLORset;

int COLORdbg_lvl(void);

int COLORgreedy (int ncount, int ecount, int *elist, int *ncolors,
        COLORset **colorclasses);

void COLORinit_set (COLORset *s);
void COLORfree_set (COLORset *s);
void COLORfree_sets (COLORset **s,int* nsets);
int  COLORcheck_set(COLORset* set, int ncount, int ecount, const int elist[]);

#define COLOR_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define COLORcheck_rval(rval,msg) {                                        \
    if ((rval)) {                                                          \
       fprintf (stderr, "%s at %s, line %d\n", (msg),__FILE__,__LINE__);   \
        goto CLEANUP;                                                      \
    }                                                                      \
}

#define COLORcheck_NULL(item,msg) {                                        \
    if ((!item)) {                                                         \
       fprintf (stderr, "%s at %s, line %d\n", (msg),__FILE__,__LINE__);   \
        rval = 1;                                                          \
        goto CLEANUP;                                                      \
    }                                                                      \
}

#endif  /* __COLOR_H */
