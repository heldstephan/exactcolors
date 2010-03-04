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
    int   rval     = 0;
    char* bosshost = (char *) NULL;
    char  task;
    char myname[MAX_PNAME_LEN],my_hostname[MAX_PNAME_LEN];
    pid_t my_pid   = getpid();
    COLORproblem  colorproblem;
    colordata* root_cd = &(colorproblem.root_cd);
    COLOR_SFILE *s = (COLOR_SFILE *) NULL;
    double cputime = COLORcpu_time();
    COLORset_dbg_lvl(-1);
    COLORproblem_init(&colorproblem);

    if (ac != 2) {
        usage (av[0]);
        rval = 1;  goto CLEANUP;
    }

    bosshost = av[1];

    rval = gethostname (my_hostname, MAX_PNAME_LEN - 1);
    COLORcheck_rval (rval, "gethostname failed");

    printf ("Machine Name: %s pid: %lld\n", my_hostname,(long long) my_pid);
    fflush (stdout);

    snprintf(myname,MAX_PNAME_LEN, "%s:%lld",
             my_hostname, (long long) my_pid);

    while (1) {
       int include_bestcolors = 0;
              
       rval = open_connection(&s, bosshost);
       COLORcheck_rval(rval,"open_connection failed.");

       rval = COLORsafe_swrite_string (s, myname);
       COLORcheck_rval (rval, "COLORsafe_swrite_string failed (NAME)");

       rval = COLORsafe_swrite_char (s, COLOR_BOSS_SEND);
       COLORcheck_rval (rval, "COLORsafe_swrite_char failed (SEND)");

       rval = COLORsafe_sread_char (s, &task);
       COLORcheck_rval (rval, "COLORsafe_sread_char failed (task)");

       if (task == COLOR_BOSS_NO) {
          COLORsafe_sclose (s);
          printf ("No more work - waiting\n");
          sleep(15);
          continue;
       }
       if (task == COLOR_BOSS_EXIT) {
          COLORsafe_sclose (s);
          printf ("No more work - shut down\n");
          goto CLEANUP;
       }

       if (task == COLOR_BOSS_YES) {
          int adopt_id = 1;
          include_bestcolors = 0;
          

          rval = receive_colordata (s, root_cd,adopt_id,include_bestcolors,
                                    &(colorproblem));
          COLORcheck_rval (rval, "receive_prob failed");

          COLORsafe_sclose (s);
       
       
          
          assert(task == COLOR_BOSS_YES);
          
          /*        sleep(3); */
          cputime = -COLORcpu_time();          
          rval = build_lp(root_cd);
          COLORcheck_rval(rval,"Failed in build_lp");
          
          rval = create_branches(root_cd,&(colorproblem));
          COLORcheck_rval(rval,"Failed in create_branches");
          cputime += COLORcpu_time();          
          
          rval = open_connection(&s, bosshost);
          COLORcheck_rval(rval,"open_connection failed.");
          
          printf ("Send the completed problem %d to boss\n", root_cd->id);
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

       }
       COLORsafe_sclose (s);
       COLORproblem_free(&colorproblem);
       COLORproblem_init(&colorproblem);
   
    }

 CLEANUP:
    COLORproblem_free(&colorproblem);
    
    return rval;
}

static void usage (char *f)
{
   fprintf (stderr, "Usage %s: bossname\n", f);
}
