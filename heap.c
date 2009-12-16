#include <stdlib.h>
#include <assert.h>

#include "heap.h"

/* #define HEAP_INTEGRITY_CHECKS 1 */

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


COLOR_MAYBE_UNUSED
static int COLORNWTheap_integrity(COLORNWTHeap* heap)
{
   int rval = 0;
   int i;
   int* perm = heap->perm;
   int* iperm = heap->iperm;

   for (i = 0 ; i < heap->end; ++i) {
      rval = ! (iperm[perm[i]] == i);
      COLORcheck_rval(rval,"Failed: iperm[perm[i]] == i");

      rval = ! (perm[iperm[i]] == i);
      COLORcheck_rval(rval,"Failed: perm[iperm[i]] == i");
      if (i > 0) {
         int parent = i >> 1;
         rval = (heap->elms[perm[parent]].key > heap->elms[perm[i]].key);
         COLORcheck_rval(rval,"Failed: heap order");
      }
   }
 CLEANUP:
   return rval;
}


#ifdef HEAP_INTEGRITY_CHECKS
#define HEAP_INTEGRITY(rval,heap,msg) {              \
      rval = COLORNWTheap_integrity(heap);           \
      COLORcheck_rval(rval,msg);                     \
}
#else
#define HEAP_INTEGRITY(rval,heap,msg)
#endif

int COLORNWTheap_init(COLORNWTHeap** heap,
                      int size)
{
   int rval = 0;
   if (size == 0) {rval = 1; goto CLEANUP;}


   *heap = (COLORNWTHeap*) COLOR_SAFE_MALLOC (1,COLORNWTHeap);
   COLORcheck_NULL(*heap,"Failed to allocate heap");


   (*heap)->perm  = (int*) NULL;
   (*heap)->iperm = (int*) NULL;
   (*heap)->elms  = (COLORNWTHeapElm*) NULL;

   (*heap)->end = 0;
   size += 2;
   (*heap)->size = size;


   (*heap)->perm = (int*) COLOR_SAFE_MALLOC(size,int);
   COLORcheck_NULL((*heap)->perm,"Failed to allocate (*heap)->perm");


   (*heap)->iperm = (int*) COLOR_SAFE_MALLOC(size,int);
   COLORcheck_NULL((*heap)->iperm,"Failed to allocate (*heap)->iperm");

   (*heap)->elms = 
      (COLORNWTHeapElm*) COLOR_SAFE_MALLOC(size,COLORNWTHeapElm);
   COLORcheck_NULL((*heap)->elms,"Failed to allocate (*heap)->elms");

   
   /* Use sentenials at beginning and end.*/
   (*heap)->elms[0].key       =   COLORNWT_MIN;
   (*heap)->perm[0] = (*heap)->iperm[0] = 0;
   
   (*heap)->elms[size-1].key  =  COLORNWT_MAX;
   (*heap)->perm[size-1] = (*heap)->iperm[size-1] = size-1;
   
   COLORNWTheap_reset(*heap);

   HEAP_INTEGRITY(rval,*heap,
                  "COLORNWTheap_integrity failed in COLORNWTheap_relabel.");
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
   assert(!COLORNWTheap_integrity(heap));
}

static int COLORNWTheap_liftup(COLORNWTHeap* heap,
                                int           pos)
{
   int swaps   = 0;
   int href     = heap->perm[pos];
   int* perm   = heap->perm;
   int* iperm  = heap->iperm;
   int  parent = pos >> 1;
   COLORNWT key = heap->elms[href].key;
   /* The sentinel at index 0 will stop the loop.*/
   while (heap->elms[perm[parent]].key > key ) {
      /* Move the parent down .*/
      perm[pos] = perm[parent];
      iperm[perm[pos]] = pos;
      pos     = parent;
      parent  >>= 1;
      ++swaps;
   }

   /* If elm at href was lifted up (to pos), update perm arrays.*/
   if(href != heap->perm[pos]) {
      perm[pos]         = href;
      iperm[href]        = pos;
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
   int  ref = perm[pos];
   int  key = heap->elms[ref].key;

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

   perm[pos] = ref;
   iperm[ref] = pos;

   return swaps;
}


int COLORNWTheap_insert (COLORNWTHeap* heap,
                         int*          href,
                         COLORNWT      key,
                         void*         obj)
{
   int rval = 0;

   (heap->end)++;
   
   if (heap->end  >= heap->size) {
      int i;

      heap->size = heap->size * 2 + 2;
      /* realloc memory */
      heap->perm = realloc(heap->perm,heap->size);
      COLORcheck_NULL(heap->perm,"Failed to reallocate heap->perm");

      heap->iperm = realloc(heap->iperm,heap->size);
      COLORcheck_NULL(heap->iperm,"Failed to reallocate heap->iperm");

      heap->elms = realloc(heap->elms,heap->size);
      COLORcheck_NULL(heap->elms,"Failed to reallocate heap->elms");
      
      for (i = heap->end; i < heap->size; ++i) {
         heap->elms[i].key = COLORNWT_MAX;
         heap->perm[i] = heap->iperm[i] = i;
      }
      assert(!COLORNWTheap_integrity(heap));
   }
   
   heap->elms[heap->perm[heap->end]].obj  = obj;
   heap->elms[heap->perm[heap->end]].key  = key;
   *href = heap->end;
   
   COLORNWTheap_liftup(heap,heap->end);

   HEAP_INTEGRITY(rval,heap,"COLORNWTheap_integrity failed in COLORNWTheap_insert.");

 CLEANUP:
   return rval;
}

int COLORNWTheap_remove (COLORNWTHeap* heap,
                         int           href)
{
   int rval = 0;
   int heap_pos = heap->iperm[href];
   
   if (heap_pos > heap->end) {
      fprintf(stderr,"COLORNWTheap_remove error: href does not exist!\n");
      rval = 1; goto CLEANUP;
   }
   
   heap->elms[href].key = COLORNWT_MAX;
   
   /* A successive uplift of the minimum child might be faster but for
      the time beeing I'm too lazy to implement this.  Instead I lift
      the last element to the hole and sift it down.
   */
   if (heap_pos < heap->end) {
      /* swap last elm to heap_pos  */
      heap->perm[heap_pos] = heap->perm[heap->end];
      heap->perm[heap->end] = href;
   
      /* swap iperm */
      heap->iperm[heap->perm[heap_pos]] = heap_pos;
      heap->iperm[heap->perm[heap->end]] = heap->end;

      
      /* Move down elm at index heap_pos. 
         It cannot travel to heap->end, as that element has
          key of COLORNWT_MAX now.
       */
      rval = COLORNWTheap_relabel(heap, 
                                  heap->perm[heap_pos], 
                                  heap->elms[heap->perm[heap_pos]].key);
   } 
   (heap->end)--;

   HEAP_INTEGRITY(rval,heap,"COLORNWTheap_integrity failed in COLORNWTheap_remove.");

 CLEANUP:
   return rval;
}

void* COLORNWTheap_min(COLORNWTHeap* heap)
{
   int href;
   void* obj;

   assert(!COLORNWTheap_integrity(heap));

   if(COLORNWTheap_empty(heap)) {
      return (void*) NULL;
   }
   
   href = heap->perm[1];
   obj = heap->elms[href].obj;

   /* swap last elm to index 1 */
   heap->perm[1] = heap->perm[heap->end];
   heap->perm[heap->end] = href;
   
   /* swap iperm */
   heap->iperm[heap->perm[1]] = 1;
   heap->iperm[heap->perm[heap->end]] = heap->end;

   
   heap->elms[heap->perm[heap->end]].key = COLORNWT_MAX;
   
   (heap->end)--;

   /* Move down elm at index 1. */
   COLORNWTheap_siftdown(heap, 1);

   assert(!COLORNWTheap_integrity(heap));

   return obj;
}

int COLORNWTheap_get_key (const COLORNWTHeap* heap,
                          int           href)
{
   assert(href < heap->size);
   return heap->elms[href].key;
}

int COLORNWTheap_size (const COLORNWTHeap* heap)
{
   return heap->end;
}




int COLORNWTheap_decrease_key (COLORNWTHeap* heap,
                               int           href,
                               COLORNWT      new_key)
{
   int rval = 0;
   int heap_pos = heap->iperm[href];

   if (heap->elms[href].key < new_key) { 
      fprintf(stderr,"COLORNWTheap_decrease_key error: new_key is greater than old key!\n");
      rval = 1; goto CLEANUP;
   }
   heap->elms[href].key = new_key;

   COLORNWTheap_liftup(heap,heap_pos);

   rval = COLORNWTheap_integrity(heap);
   COLORcheck_rval(rval,"COLORNWTheap_integrity failed in COLORNWTheap_decrease_key.\n");

 CLEANUP:
   return rval;
}

int COLORNWTheap_relabel (COLORNWTHeap* heap,
                          int           href,
                          COLORNWT      new_key)
{
   int rval = 0;
   int heap_pos = heap->iperm[href];
   heap->elms[href].key = new_key;
   if (!COLORNWTheap_liftup(heap,heap_pos)) {
      COLORNWTheap_siftdown(heap,heap_pos);
   }

   HEAP_INTEGRITY(rval,heap,"COLORNWTheap_integrity failed in COLORNWTheap_relabel.");

#ifdef HEAP_INTEGRITY_CHECKS   
 CLEANUP:
#endif
   return rval;
}
