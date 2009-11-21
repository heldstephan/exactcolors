#ifndef __COLOR_H
#define __COLOR_H

typedef struct COLORset {
    int count;
    int *members;
    struct COLORset *next;
} COLORset;

int COLORdbg_lvl(void);

int COLORgreedy (int ncount, int ecount, int *elist, int *ncolors,
        COLORset **colorclasses);

void *COLORutil_allocrus (size_t size);
void COLORutil_freerus (void *p);
void COLORinit_set (COLORset *s);
void COLORfree_set (COLORset *s);
void COLORfree_sets (COLORset **s,int* nsets);
int  COLORcheck_set(COLORset* set, int ncount, int ecount, const int elist[]);

double COLORwall_time (void);
double COLORcpu_time (void);


#define COLOR_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define COLOR_SAFE_MALLOC(nnum,type)                                       \
    (type *) COLORutil_allocrus (((size_t) (nnum)) * sizeof (type))

#define COLOR_FREE(object,type) {                                          \
    COLORutil_freerus ((void *) (object));                                 \
    object = (type *) NULL;                                                \
}

#define COLOR_IFFREE(object,type) {                                        \
    if ((object)) COLOR_FREE ((object),type);                              \
}

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
