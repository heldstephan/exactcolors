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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <sys/resource.h>
#include <sys/stat.h>

#include "color_defs.h"
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

int COLORfile_exists(const char* filename)
{
   FILE * file = fopen(filename, "r");
   if (file)
   {
      fclose(file);
      return 1;
   }
   return 0;
}

int COLORdir_exists(const char* dirname)
{
   struct stat st;

   return (stat(dirname,&st) == 0);
}


int COLORdir_create(const char* dirname)
{
   int    prval = 0;
   int    rval  = 0;

   prval = mkdir(dirname,(S_IRUSR | S_IWUSR | S_IXUSR));
   COLORcheck_fileio(prval,"Failed to mkdir");
   
 CLEANUP:
   
   return rval;
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
    for (i = 0; i < 165; i++) COLORutil_lprand (r);
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
        s->age   = 0;
        s->next  = (COLORset*) NULL;
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
      for (i = 0; i < *nsets; i++) {
         COLORfree_set (& (*sets)[i]);
      }
      COLOR_FREE (*sets, COLORset);
   }
   *nsets = 0;
}

int COLORcopy_sets (COLORset **s,int *nsets,
                    const COLORset *src_s, int src_nsets)
{
   int rval = 0;
   int i;
   *nsets = src_nsets;
   *s = (COLORset*) COLOR_SAFE_MALLOC(src_nsets,COLORset);
   COLORcheck_NULL(*s,"Failed to allocate *s");

   for (i = 0; i < src_nsets; ++i) {
      (*s)[i].count = src_s[i].count;
      (*s)[i].members = (int*) COLOR_SAFE_MALLOC(src_s[i].count,int);
      COLORcheck_NULL((*s)[i].members,"Failed to allocate (*s)[i].members");
      memcpy((*s)[i].members,src_s[i].members,src_s[i].count * sizeof(int));
   }

 CLEANUP:
   return rval;
}


double COLORwall_time (void)
{
    return (double) time (0);
}

double COLORcpu_time (void)
{
    struct rusage ru;
    double t;

    getrusage (RUSAGE_SELF, &ru);

    t = ((double) ru.ru_utime.tv_sec) +
        ((double) ru.ru_utime.tv_usec) / 1000000.0;
    return t;
}
