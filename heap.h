#ifndef __COLOR_HEAP_H
#define __COLOR_HEAP_H

#include "color_defs.h"


typedef struct COLORNWTHeap_t COLORNWTHeap;

int COLORNWTheap_init(COLORNWTHeap** heap,
                      int            size);

int COLORNWTheap_free(COLORNWTHeap*  heap);

void COLORNWTheap_reset(COLORNWTHeap* heap);


int COLORNWTheap_insert (COLORNWTHeap* heap,
                         int*          href,
                         COLORNWT      key,
                         void*         obj);

int COLORNWTheap_remove (COLORNWTHeap* heap,
                         int           href);


void* COLORNWTheap_min (COLORNWTHeap* heap);


int COLORNWTheap_get_key (const COLORNWTHeap* heap,
                          int           href);

int COLORNWTheap_size (const COLORNWTHeap* heap);


int COLORNWTheap_decrease_key (COLORNWTHeap* heap,
                               int           href,
                               COLORNWT      new_key);


int COLORNWTheap_relabel (COLORNWTHeap* heap,
                          int           href,
                          COLORNWT      new_key);



#endif
