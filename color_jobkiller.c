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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>
#include <assert.h>
#include <limits.h>

#include "color_defs.h"
#include "graph.h"
#include "color.h"
#include "color_private.h"

#include "bbsafe.h"

int main (int ac, char **av);
static void usage (char *f);

static int open_connection(COLOR_SFILE** s,const char* bosshost)
{
   int rval = 0;
   int k = 0;
   do {
      *s = COLORsafe_snet_open (bosshost, COLOR_BOSS_PORT);
      if (! (*s)) {
         fprintf (stderr, "COLORsafe_snet_open failed\n");
         sleep (10);
      }
      k++;
   } while (! (*s) && k < 10000);

   if (! (*s)) {
      fprintf (stderr, "Could not connect in %d trys.\n", k);
      goto CLEANUP;
   }
 CLEANUP:
   return rval;
}


int main (int ac, char **av)
{
    int   rval     = 0;
    char myname[MAX_PNAME_LEN],my_hostname[MAX_PNAME_LEN];
    pid_t my_pid   = getpid();
    int  cd_id;
    COLOR_SFILE *s = (COLOR_SFILE *) NULL;

    rval = COLORprogram_header (ac,av);
    COLORcheck_rval(rval, "Failed in COLORprogram_header");
    
    
    if (ac != 2) {
        usage (av[0]);
        rval = 1;  goto CLEANUP;
    }


    cd_id = atoi(av[1]);

    rval = gethostname (my_hostname, MAX_PNAME_LEN - 1);
    COLORcheck_rval (rval, "gethostname failed");

    snprintf(myname,MAX_PNAME_LEN, "%s:%lld",
             my_hostname, (long long) my_pid);

    printf("Trying to remove job id = %d from color job on %s.\n",
           cd_id, my_hostname);
      
    rval = open_connection(&s, my_hostname);
    COLORcheck_rval(rval,"open_connection failed.");
       
    rval = COLORsafe_swrite_string (s, myname);
    COLORcheck_rval (rval, "COLORsafe_swrite_string failed (NAME)");

    rval = COLORsafe_swrite_char (s, COLOR_BOSS_REMOVE_JOB);
    COLORcheck_rval (rval, "COLORsafe_swrite_char failed (REMOVE_JOB)");

    rval = COLORsafe_swrite_int (s, cd_id);
    COLORcheck_rval (rval, "COLORsafe_swrite_int failed (cd_id)");

    COLORsafe_sclose (s);

 CLEANUP:

    return rval;
}

static void usage (char *f)
{
   fprintf (stderr, "Usage %s: bossname\n", f);
}
