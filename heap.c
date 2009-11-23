#include <stdlib.h>
#include <assert.h>

#include "heap.h"

typedef struct COLORNWTHeapElm {
   COLORNWT key;
   void*    obj;
} COLORNWTHeapElm;

struct COLORNWTHeap_t {
   int     end;
   int     size;
   int*    perm;
   int*    iperm;
   
   /* elements at index 0 and size-1 are sentinels.
      The actual elements are stored at indices 1..end < size -1.
   */
   COLORNWTHeapElm* elms;
};

static int COLORNWTheap_empty(COLORNWTHeap* heap)
{
   assert(heap->end >= 0);
   return !(heap->end);
}

int COLORNWTheap_init(COLORNWTHeap** heap,
                      int size)
{
   int rval = 0;
   if (size == 0) {rval = 1; goto CLEANUP;}


   *heap = (COLORNWTHeap*) malloc(sizeof(COLORNWTHeap));
   COLORcheck_NULL(*heap,"Failed to allocate heap");


   (*heap)->perm  = (int*) NULL;
   (*heap)->iperm = (int*) NULL;
   (*heap)->elms  = (COLORNWTHeapElm*) NULL;

   (*heap)->end = 0;
   size += 2;
   (*heap)->size = size;


   (*heap)->perm = (int*) malloc(size * sizeof(int));
   COLORcheck_NULL((*heap)->perm,"Failed to allocate (*heap)->perm");


   (*heap)->iperm = (int*) malloc(size * sizeof(int));
   COLORcheck_NULL((*heap)->iperm,"Failed to allocate (*heap)->iperm");

   (*heap)->elms = 
      (COLORNWTHeapElm*) malloc(size * sizeof(COLORNWTHeapElm));
   COLORcheck_NULL((*heap)->elms,"Failed to allocate (*heap)->elms");

   
   /* Use sentenials at beginning and end.*/
   (*heap)->elms[0].key       =   COLORNWT_MIN;
   (*heap)->perm[0] = (*heap)->iperm[0] = 0;
   
   (*heap)->elms[size-1].key  =  COLORNWT_MAX;
   (*heap)->perm[size-1] = (*heap)->iperm[size-1] = size-1;
   
   COLORNWTheap_reset(*heap);
 CLEANUP:
   if(rval) {
      COLORNWTheap_free(*heap);
      *heap = (COLORNWTHeap*) NULL;
   }
   return rval;
}

int COLORNWTheap_free(COLORNWTHeap*  heap)
{
   if (heap) { 
      if (heap->perm) free(heap->perm);
      if (heap->iperm) free(heap->iperm);
      if (heap->elms) free(heap->elms);
      free (heap);
      
   }
   return 0;
}

void COLORNWTheap_reset(COLORNWTHeap* heap)
{
   int i;

   heap->end = 0;

   for(i = 1; i + 1 < heap->size; ++i){
      heap->elms[i].key = COLORNWT_MAX;
      heap->perm[i] = heap->iperm[i] = i;
   }
}

static int COLORNWTheap_liftup(COLORNWTHeap* heap,
                                int           pos)
{
   int swaps  = 0;
   int idx    = heap->perm[pos];
   int* perm  = heap->perm;
   int* iperm = heap->iperm;
   int pred   = pos >> 1;
   COLORNWT key = heap->elms[pos].key;

   /* The sentinel at index 0 will stop the loop.*/
   while (heap->elms[perm[pred]].key > key ) {
      /* Move the parent down .*/
      perm[pos] = perm[pred];
      iperm[perm[pos]] = pos;
      pos     = pred;
      pred  >>= 1;
      ++swaps;
   }

   /* If elm at idx was lifted up (to pos), update perm arrays.*/
   if(idx != heap->perm[pos]) {
      perm[pos]         = idx;
      iperm[idx]        = pos;
   }

   return swaps;
}

static int COLORNWTheap_siftdown(COLORNWTHeap* heap,
                                 int           pos)
{
   int swaps    = 0;
   int end_half = heap->end / 2;

   int* perm  = heap->perm;
   int* iperm = heap->iperm;
   COLORNWTHeapElm* elms = heap->elms;
   int  minc,rightc;
   int  idx = perm[pos];
   int  key = heap->elms[idx].key;

   while (pos <= end_half) {
      minc  = pos << 1;  /* j = k*2 */
      rightc = minc + 1;
      /* set minc to minimum of left and right child */
      if ( elms[perm[minc]].key  >  elms[perm[rightc]].key  )
         minc = rightc;
      
      if ( key <= elms[perm[minc]].key )
         break;   
      
      /* move element 'minc' up in the heap */
      perm[pos] = perm[minc];
      iperm[perm[pos]] = pos;
      ++swaps;
      pos = minc;
   }

   perm[pos] = idx;
   iperm[idx] = pos;

   return swaps;
}


int COLORNWTheap_insert (COLORNWTHeap* heap,
                         int*          pos,
                         COLORNWT      key,
                         void*         obj)
{
   int rval = 0;

   (heap->end)++;
   
   assert (heap->end  < heap->size);
   
   heap->elms[heap->perm[heap->end]].obj  = obj;
   heap->elms[heap->perm[heap->end]].key  = key;
   *pos = heap->end;
   
   COLORNWTheap_liftup(heap,heap->end);

   return rval;
}

void* COLORNWTheap_min(COLORNWTHeap* heap)
{
   int idx;
   void* obj;
   if(COLORNWTheap_empty(heap)) {
      return (void*) NULL;
   }
   
   idx = heap->perm[1];
   obj = heap->elms[idx].obj;

   /* swap last elm to index 1 */
   heap->perm[1] = heap->perm[heap->end];
   heap->perm[heap->end] = idx;
   
   /* swap iperm */
   heap->iperm[heap->perm[1]] = 1;
   heap->iperm[heap->perm[heap->end]] = heap->end;

   
   heap->elms[heap->perm[heap->end]].key = COLORNWT_MAX;
   
   (heap->end)--;

   /* Move down elm at index 1. */
   COLORNWTheap_siftdown(heap, 1);

   return obj;
}

int COLORNWTheap_decrease_key (COLORNWTHeap* heap,
                               int           pos,
                               COLORNWT      new_key)
{
   heap->elms[pos].key = new_key;

   COLORNWTheap_liftup(heap,heap->iperm[pos]);
   return 0;
}

