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
#include "color.h"

#include "plotting.h"

typedef unsigned char RGBColor_t[3];

#define COLORS 16
const char* ColorTable[COLORS] = 
   {        "red",   "yellow",          "blue",        "brown",
         "orange",    "green",          "cyan",      "magenta",
     "darksalmon",    "gold1", "darkslateblue",   "darkviolet",
     "darkorange","darkgreen", "darkturquoise", "midnightblue"
   };


int COLORplot_graphviz(const char* filename,
                       int ncount, int ecount, const int elist[], 
                       int sets[])
{
   int rval = 0,i;


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
