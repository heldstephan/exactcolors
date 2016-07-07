#ifndef __COLOR_DEFS_H
#define __COLOR_DEFS_H
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

#include <stdio.h>
#include <limits.h>
#ifndef COMPILE_FOR_VALGRIND
#include "rounding_mode.h"
#else 
#include <float.h>
#include <math.h>
#endif

#ifdef __GNUC__
    #define COLOR_MAYBE_UNUSED __attribute__((used))
#else
    #define COLOR_MAYBE_UNUSED
#endif


#define COLOR_MAXINT (2147483647)
#define MAX_PNAME_LEN 128


double COLORwall_time (void);
double COLORcpu_time (void);

void *COLORutil_allocrus (size_t size);
void  COLORutil_freerus (void *p);

int COLORfile_exists(const char* filename);
int COLORdir_exists(const char* dirname);
int COLORdir_create(const char* dirname);

int COLORdbg_lvl(void);
void COLORset_dbg_lvl(int dbglvl);


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
       fflush(stdout);                                                     \
       fprintf (stderr, "%s at %s, line %d\n", (msg),__FILE__,__LINE__);   \
       goto CLEANUP;                                                       \
    }                                                                      \
}

#define COLORcheck_fileio(rval,msg) {                                      \
    if (prval < 0) {                                                       \
       fflush(stdout);                                                     \
       fprintf (stderr, "%s at %s, line %d\n", (msg),__FILE__,__LINE__);   \
       rval = 1; goto CLEANUP;                   			   \
    }                                                                      \
}

#define COLORcheck_NULL(item,msg) {                                        \
    if ((!item)) {                                                         \
       fflush(stdout);                                                     \
       fprintf (stderr, "%s at %s, line %d\n", (msg),__FILE__,__LINE__);   \
       rval = 1;                                                           \
       goto CLEANUP;                                                       \
    }                                                                      \
}



/* The type of the node weights in upcoming
   maximum weighted independent set and clique problems.
   COLORNWT = COLOR Node Weight Type.
*/
typedef int COLORNWT;
#define COLORNWT_MAX INT_MAX
#define COLORNWT_MIN INT_MIN

COLOR_MAYBE_UNUSED static inline int COLORnode_comparator(const void* v1,const void* v2)
{
   int i1 = *(const int*) v1;
   int i2 = *(const int*) v2;

   return i1-i2;
}

COLOR_MAYBE_UNUSED static inline double COLORDBLmax(double a,double b)
{
      return a > b ? a : b;
}

COLOR_MAYBE_UNUSED static inline double COLORDBLmin(double a,double b)
{
      return a < b ? a : b;
}


COLOR_MAYBE_UNUSED static inline COLORNWT COLORNWTmax(COLORNWT a,COLORNWT b)
{
      return a > b ? a : b;
}

COLOR_MAYBE_UNUSED static inline COLORNWT COLORNWTmin(COLORNWT a,COLORNWT b)
{
      return a < b ? a : b;
}

#ifndef COMPILE_FOR_VALGRIND
COLOR_MAYBE_UNUSED static double COLORsafe_lower_dbl(COLORNWT numerator,COLORNWT denominator)
{
   double result;
   int dummy_val;

   int    oldround = COLOR_get_rounding(&dummy_val);
   double denom_mult;
   COLOR_set_round_up(&dummy_val);
   denom_mult = denominator;

   COLOR_set_round_down(&dummy_val);
   denom_mult = 1 / denom_mult;

   result = (double) numerator * denom_mult;
   COLOR_set_rounding(oldround,&dummy_val);

   return result;
}
#else
COLOR_MAYBE_UNUSED static double COLORsafe_lower_dbl(COLORNWT numerator,COLORNWT denominator)
{
   double result;
   double denom_mult;
   denom_mult = denominator;
   denom_mult = nextafter(denom_mult, DBL_MAX);

   denom_mult = 1 / denom_mult;
   denom_mult = nextafter(denom_mult, -DBL_MAX);

   result = (double) numerator * denom_mult;
   result = nextafter(result, -DBL_MAX);

   return result;
}
#endif

COLOR_MAYBE_UNUSED static double COLORunsafe_dbl(COLORNWT numerator,COLORNWT denominator)
{
   return  (double) numerator / (double) denominator;
}



#endif
