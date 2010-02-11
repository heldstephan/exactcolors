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
#include <getopt.h>
#include <string.h>
#include <fenv.h>
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
/*          sleep (100); */
         sleep (5);
      }
      k++;
   } while (! (*s) && k < 2 /* 5 */);
        
   if (! (*s)) {
      fprintf (stderr, "Could not connect in %d trys.\n", k);
      goto CLEANUP;
   }
 CLEANUP:
   return rval;
}
                           

int main (int ac, char **av)
{
    int rval = 0;
    char* bosshost = (char *) NULL;
    char task, myname[256];
    colordata   root_cdata;
    colordata*  root_cd = &root_cdata;
    COLOR_SFILE *s = (COLOR_SFILE *) NULL;
    double cputime = COLORcpu_time();
    COLORset_dbg_lvl(-1);

    if (ac != 2) {
        usage (av[0]);
        rval = 1;  goto CLEANUP;
    }

    bosshost = av[1];

    rval = gethostname (myname, 127);
    COLORcheck_rval (rval, "gethostname failed");
    printf ("Machine Name: %s\n", myname);  fflush (stdout);

    while (1) {
       int include_bestcolors = 0;
       
       init_colordata(root_cd);
       
       rval = open_connection(&s, bosshost);
       COLORcheck_rval(rval,"open_connection failed.");

       rval = COLORsafe_swrite_string (s, myname);
       COLORcheck_rval (rval, "COLORsafe_swrite_string failed (NAME)");

       rval = COLORsafe_swrite_char (s, COLOR_BOSS_SEND);
       COLORcheck_rval (rval, "COLORsafe_swrite_char failed (SEND)");

       rval = COLORsafe_sread_char (s, &task);
       COLORcheck_rval (rval, "COLORsafe_sread_char failed (task)");

       if (task == COLOR_BOSS_YES) {
          include_bestcolors = 0;
          int adopt_id = 1;
          rval = receive_colordata (s, root_cd,adopt_id,include_bestcolors);
          COLORcheck_rval (rval, "receive_prob failed");
       } 

       COLORsafe_sclose (s);

       if (task == COLOR_BOSS_NO) {
          printf ("No more work - waiting\n");
          sleep(15);
          continue;
       }
       if (task == COLOR_BOSS_EXIT) {
          printf ("No more work - shut down\n");
          goto CLEANUP;
       }
       

       assert(task == COLOR_BOSS_YES);

/*        sleep(3); */
       cputime = -COLORcpu_time();          
       rval = build_lp(root_cd);
       COLORcheck_rval(rval,"Failed in build_lp");

       rval = create_branches(root_cd);
       COLORcheck_rval(rval,"Failed in create_branches");
       cputime += COLORcpu_time();          

       rval = open_connection(&s, bosshost);
       COLORcheck_rval(rval,"open_connection failed.");

       printf ("Send boss the completed problem %d\n", root_cd->id);
       fflush (stdout);

       rval = COLORsafe_swrite_string (s, myname);
       COLORcheck_rval (rval, "COLORsafe_swrite_string failed (NAME)");

       rval = COLORsafe_swrite_char (s, COLOR_BOSS_RECEIVE);
       COLORcheck_rval (rval, "COLORsafe_swrite_char failed (RECEIVE)");

       rval = COLORsafe_swrite_int (s, root_cd->id);
       COLORcheck_rval (rval, "COLORsafe_swrite_int failed (root_cd->id)");

       rval = COLORsafe_swrite_double (s, cputime);
       COLORcheck_rval (rval, "COLORsafe_swrite_double failed (cputime)");

       include_bestcolors = 1;
       rval = send_colordata (s, root_cd,include_bestcolors);
       COLORcheck_rval (rval, "send_result failed");

       COLORsafe_sclose (s);

       free_colordata(root_cd);

    }


 CLEANUP:

    return rval;
}

static void usage (char *f)
{
   fprintf (stderr, "Usage %s: bossname\n", f);
}
