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
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int node_id(int i, int j, int m)
{
   return m * (i - 1) + j;
}

static void print_node_edges(int i, int j, int n, int m)
{
   int k, l;
   
   int v = node_id(i,j,m);

   /* go north-east */
   for(k = 1; (i - k > 0) && (j + k <= m); ++k) {
      int w = node_id(i - k, j + k, m);
      printf("e %d %d\n", v, w);
   }
   

   /* go east */
   for(l = j + 1; l <= m; ++l) {
      int w = node_id(i,l,m);
      printf("e %d %d\n", v, w);
   }

   /* go south-east */
   for(k = 1; (i + k <= n) && (j + k <= m); ++k) {
      int w = node_id(i + k,j + k,m);
      printf("e %d %d\n", v, w);
   }


   /* go south */
   for(k = i + 1; k <= n; ++k) {
      int w = node_id(k, j,m);
      printf("e %d %d\n", v, w);
   }
}

int main(int argc, char** args)
{
   int n,m;
   int i,j;
   int ecount = 0;
   int ncount = 0;
   int pcount = 0;
   int dcount = 0;
   int mindim = 0;
   int maxdim = 0;
   n = m = 5;


   if (argc == 1) {
      printf("Usage: queen <number of rows> <number of columns>\n");
      printf("If the number of columns is omitted, the number of columns will be set to the number of rows.\n");
      exit (0);
   }
   if (argc >=2) {
      n = atoi(args[1]);
   }
   if (argc >=3) {
      m = atoi(args[2]);
   } else {
      m = n;
   }
   ncount = n*m;
   mindim = (n < m) ? n : m;
   maxdim = (n < m) ? m : n;
   pcount = m * n * (n - 1) / 2 + n * m * (m-1) / 2 ; /* axis-parallel edges */

   for (i = 1; i < mindim; i++) {
      dcount += i * (i-1) * 2;
   }
   dcount +=  (maxdim - mindim + 1) * mindim * (mindim - 1);   /* diagonal edges */
   

   ecount = dcount + pcount;
/*    printf("%d = %d + %d\n",ecount, pcount,dcount); */
      
   printf("p edges %d %d\n",ncount, ecount);

   for (i = 1; i <= n; ++i) {
      for (j = 1; j <= m; ++j) {
         print_node_edges(i,j,n,m);
      }
   }
   
   if (argc > 1) {
      n = atoi(args[1]);
   }
}
