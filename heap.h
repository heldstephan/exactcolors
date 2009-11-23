#ifndef __COLOR_HEAP_H
#define __COLOR_HEAP_H

#include "color_defs.h"


typedef struct COLORNWTHeap_t COLORNWTHeap;

int COLORNWTheap_init(COLORNWTHeap** heap,
                      int size);

int COLORNWTheap_free(COLORNWTHeap*  heap);

void COLORNWTheap_reset(COLORNWTHeap* heap);


int COLORNWTheap_insert (COLORNWTHeap* heap,
                         int*          pos,
                         COLORNWT      key,
                         void*         obj);


void* COLORNWTheap_min (COLORNWTHeap* heap);




int COLORNWTheap_decrease_key (COLORNWTHeap* heap,
                               int           pos,
                               COLORNWT      new_key);



#endif
