#include <stdio.h>

#include "color.h"

typedef unsigned char RGBColor_t[3];

int COLORplot_graphviz(int ncount, int ecount, const int elist[], int sets[]);


#define COLORS 16
const char* ColorTable[COLORS] = 
   {        "red",   "yellow",          "blue",        "brown",
         "orange",    "green",          "cyan",      "magenta",
     "darksalmon",    "gold1", "darkslateblue",   "darkviolet",
     "darkorange","darkgreen", "darkturquoise", "midnightblue"
   };


int COLORplot_graphviz(int ncount, int ecount, const int elist[], int sets[])
{
   int rval = 0,i;

   const char* filename = "graph.dot";

   FILE* file = fopen(filename,"w");
   if (!file) {
      fprintf(stderr,"Couldn't open %s for writing.\n",filename);
      rval = 1; goto CLEANUP;
   }

   fprintf(file,"graph G {\n");
   for (i = 0; i < ecount;++i) {
      fprintf(file,"    %d -- %d;\n",elist[2*i],elist[2*i+1]);
   }


   for (i = 0; i < ncount;++i) {
      if (sets) {
         int c = sets[i] % COLORS;
         printf("node %d color %d;\n",i, c);
         fprintf(file,"    %d [style=filled,fillcolor=%s];\n",
                 i,ColorTable[c]);
      }
   }

   
   fprintf(file,"}\n");
   

 CLEANUP:
   if (file) {
      fclose(file);
   }
   return rval;
}
