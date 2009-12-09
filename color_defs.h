#ifndef __COLOR_BASIS_H
#define __COLOR_BASIS_H

#include <stdio.h>
#include <limits.h>
#include <fenv.h>

#ifdef __GNUC__
    #define MAYBE_UNUSED __attribute__((used))
#else
    #define MAYBE_UNUSED
#endif


#define COLOR_MAXINT (2147483647)

double COLORwall_time (void);
double COLORcpu_time (void);

void *COLORutil_allocrus (size_t size);
void COLORutil_freerus (void *p);

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
        goto CLEANUP;                                                      \
    }                                                                      \
}

#define COLORcheck_NULL(item,msg) {                                        \
    if ((!item)) {                                                         \
       fflush(stdout);                                                     \
       fprintf (stderr, "%s at %s, line %d\n", (msg),__FILE__,__LINE__);   \
        rval = 1;                                                          \
        goto CLEANUP;                                                      \
    }                                                                      \
}



/* The type of the node weights in upcoming
   maximum weighted independent set and clique problems.
   COLORNWT = COLOR Node Weight Type.
*/
typedef int COLORNWT;
#define COLORNWT_MAX INT_MAX
#define COLORNWT_MIN INT_MIN

MAYBE_UNUSED static inline double dblmax(double a,double b) 
{
      return a > b ? a : b;
}

MAYBE_UNUSED static inline double dblmin(double a,double b) 
{
      return a < b ? a : b;
}


MAYBE_UNUSED static inline COLORNWT COLORNWTmax(COLORNWT a,COLORNWT b) 
{
      return a > b ? a : b;
}

MAYBE_UNUSED static inline COLORNWT COLORNWTmin(COLORNWT a,COLORNWT b) 
{
      return a < b ? a : b;
}

MAYBE_UNUSED static double COLORsafe_lower_dbl(COLORNWT numerator,COLORNWT denominator)
{
   double result;
   int oldround = fegetround();

   fesetround(FE_UPWARD);
   double denom_mult  = denominator;

   fesetround(FE_DOWNWARD);
   denom_mult = 1 / denom_mult;

   result = (double) numerator * denom_mult;
   fesetround(oldround);
   return result;
}

MAYBE_UNUSED static double COLORunsafe_dbl(COLORNWT numerator,COLORNWT denominator)
{
   return  (double) numerator / (double) denominator;
}


#endif
