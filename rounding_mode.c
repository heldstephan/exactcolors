#include "rounding_mode.h"

#ifndef COMPILE_FOR_VALGRIND
#include <fenv.h>
#define UNUSED(x) (void)(x)


void COLOR_set_round_up(int* dummy_val)
{
   UNUSED(dummy_val);
   fesetround(FE_UPWARD);
}

void COLOR_set_round_down(int* dummy_val)
{
   UNUSED(dummy_val);
   fesetround(FE_DOWNWARD);
}

void COLOR_set_rounding(int rounding,int* dummy_val)
{
   UNUSED(dummy_val);
   fesetround(rounding);
}
int COLOR_get_rounding(int* dummy_val)
{
   UNUSED(dummy_val);
   return fegetround();
}
#endif
