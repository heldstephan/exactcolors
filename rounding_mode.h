#ifndef ROUNDING_MODE
#define ROUNDING_MODE

/**
   The following are wrappers around fesetround.
   The dummy_val variable has to be passed, to make sure 
   that the function is called where intended and not moved
   to other places by the optimizer.
 */

void COLOR_set_round_up(int* dummy_val);
void COLOR_set_round_down(int* dummy_val);
void COLOR_set_rounding(int rounding, int* dummy_val);

int COLOR_get_rounding(int* dummy_val);

#endif
