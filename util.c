#include <stdlib.h>
#include <stdio.h>
#include <sys/resource.h>

#include "color.h"

void *COLORutil_allocrus (size_t size)
{
    void *mem = (void *) NULL;

    if (size == 0) {
        fprintf (stderr, "Warning: 0 bytes allocated\n");
    }

    mem = (void *) malloc (size);
    if (mem == (void *) NULL) {
        fprintf (stderr, "Out of memory. Asked for %d bytes\n", (int) size);
    }
    return mem;
}

void COLORutil_freerus (void *p)
{
    if (!p) {
        fprintf (stderr, "Warning: null pointer freed\n");
        return;
    }
    free (p);
}


/****************************************************************************/
/*    RNG Based on DIMACS Code.                                             */
/*                                                                          */
/*    DIMACS NOTES: This is a COMPLETELY PORTABLE generator. It will        */
/*    give IDENTICAL sequences of random numbers for any architecture with  */
/*    at least 30-bit integers, regardless of the integer representation,   */
/*    INT_MAX value, or roundoff/truncation method, etc.                    */
/*        This Truly Remarkable RNG is described more fully in              */
/*    J. Bentley's column, ``The Software Exploratorium ''. It is based on  */
/*    one in Knuth, Vol 2, Section 3.2.2 (Algorithm A).                     */
/****************************************************************************/


void COLORutil_sprand (int seed, COLORrandstate *r)
{
    int i, ii;
    int last, next;
    int *arr = r->arr;

    seed %= COLOR_PRANDMAX;
    if (seed < 0) seed += COLOR_PRANDMAX;

    arr[0] = last = seed;
    next = 1;
    for (i = 1; i < 55; i++) {
        ii = (21 * i) % 55;
        arr[ii] = next;
        next = last - next;
        if (next < 0)
            next += COLOR_PRANDMAX;
        last = arr[ii];
    }
    r->a = 0;
    r->b = 24;
    for (i = 0; i < 165; i++)
        last = COLORutil_lprand (r);
}

int COLORutil_lprand (COLORrandstate *r)
{
    int t;

    if (r->a-- == 0) r->a = 54;
    if (r->b-- == 0) r->b = 54;

    t = r->arr[r->a] - r->arr[r->b];

    if (t < 0) t += COLOR_PRANDMAX;

    r->arr[r->a] = t;

    return t;
}

double COLORutil_zeit (void)
{
    struct rusage ru;

    getrusage (RUSAGE_SELF, &ru);

    return ((double) ru.ru_utime.tv_sec) +
           ((double) ru.ru_utime.tv_usec) / 1000000.0;
}

void COLORutil_quicksort (int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1) return;

    COLOR_SWAP (len[0], len[(n - 1)/2], temp);

    i = 0; j = n; t = len[0];

    while (1) {
        do i++; while (i < n && len[i] < t);
        do j--; while (len[j] > t);
        if (j < i) break;
        COLOR_SWAP (len[i], len[j], temp);
    }
    COLOR_SWAP (len[0], len[j], temp);

    COLORutil_quicksort (len, j);
    COLORutil_quicksort (len + i, n - i);
}

void COLORutil_quicksort_reverse (int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1) return;

    COLOR_SWAP (len[0], len[(n - 1)/2], temp);

    i = 0; j = n; t = len[0];

    while (1) {
        do i++; while (i < n && len[i] > t);
        do j--; while (len[j] < t);
        if (j < i) break;
        COLOR_SWAP (len[i], len[j], temp);
    }
    COLOR_SWAP (len[0], len[j], temp);

    COLORutil_quicksort_reverse (len, j);
    COLORutil_quicksort_reverse (len + i, n - i);
}

void COLORutil_perm_quicksort (int *perm, int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1) return;

    COLOR_SWAP (perm[0], perm[(n - 1)/2], temp);

    i = 0; j = n; t = len[perm[0]];

    while (1) {
        do i++; while (i < n && len[perm[i]] < t);
        do j--; while (len[perm[j]] > t);
        if (j < i) break;
        COLOR_SWAP (perm[i], perm[j], temp);
    }
    COLOR_SWAP (perm[0], perm[j], temp);

    COLORutil_perm_quicksort (perm, len, j);
    COLORutil_perm_quicksort (perm + i, len, n - i);
}

void COLORinit_set (COLORset *s) 
{
    if (s) {
        s->members = (int *) NULL;
        s->count = 0;
    }
}

void COLORfree_set (COLORset *s) 
{
    if (s && s->members) {
        free (s->members);
        s->members = (int *) NULL;
        s->count = 0;
    }
}

void COLORfree_sets(COLORset** sets,int* nsets)
{
   int i;
   if (*sets) {
      for (i = 0; i < *nsets; i++) COLORfree_set (& (*sets)[i]);
      free (*sets);
   }
   *sets  = (COLORset*) NULL;
   *nsets = 0;
}

