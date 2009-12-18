#ifndef __COLOR_HEAP_H
#define __COLOR_HEAP_H
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
