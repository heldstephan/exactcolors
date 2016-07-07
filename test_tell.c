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
#include "bbsafe.h"

static int debug = 0;

int main (int ac, char **av);
static void usage (char *f);

int main (int ac, char **av)
{
    int k, rval = 0;
    char myname[128], *bosshost = (char *) NULL;
    COLOR_SFILE *s = (COLOR_SFILE *) NULL;

    if (ac != 2) {
        usage (av[0]);
        rval = 1;  goto CLEANUP;
    }

    bosshost = av[1];

    rval = gethostname (myname, 127);
    COLORcheck_rval (rval, "gethostname failed");
    printf ("Machine Name: %s\n", myname);  fflush (stdout);
    printf ("Tell Boss %s to Exit\n",  bosshost);  fflush (stdout);

    k = 0;
    do {
        s = COLORsafe_snet_open (bosshost, COLOR_BOSS_PORT);
        if (!s) {
            fprintf (stderr, "COLORsafe_snet_open failed\n");
            sleep (100);
        }
        k++;
    } while (!s && k < 5);

    if (!s) {
        fprintf (stderr, "Could not connect in %d trys.\n", k);
        goto CLEANUP;
    }

    rval = COLORsafe_swrite_string (s, myname);
    COLORcheck_rval (rval, "COLORsafe_swrite_string failed (NAME)");

    rval = COLORsafe_swrite_char (s, COLOR_BOSS_EXIT);
    COLORcheck_rval (rval, "COLORsafe_swrite_char failed (EXIT)");

    COLORsafe_sclose (s);

CLEANUP:

    return rval;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage %s: bossname\n", f);
}

int COLORdbg_lvl() {
   return debug;
}
