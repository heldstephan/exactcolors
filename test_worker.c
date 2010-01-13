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
#include "bbsafe.h"

static int debug = 0;

int main (int ac, char **av);
static int receive_prob (COLOR_SFILE *s, int *probid, int *ncount, int *ecount,
        int **pelist);
static int send_result (COLOR_SFILE *s, int probid);
static void usage (char *f);

int main (int ac, char **av)
{
    int k, probid = 0, rval = 0;
    int ncount, ecount;
    int *elist = (int *) NULL;
    char task, myname[128], *bosshost = (char *) NULL;
    COLOR_SFILE *s = (COLOR_SFILE *) NULL;

    if (ac != 2) {
        usage (av[0]);
        rval = 1;  goto CLEANUP;
    }

    bosshost = av[1];

    rval = gethostname (myname, 127);
    COLORcheck_rval (rval, "gethostname failed");
    printf ("Machine Name: %s\n", myname);  fflush (stdout);

    while (1) {
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

        rval = COLORsafe_swrite_char (s, COLOR_BOSS_SEND);
        COLORcheck_rval (rval, "COLORsafe_swrite_char failed (SEND)");

        rval = COLORsafe_sread_char (s, &task);
        COLORcheck_rval (rval, "COLORsafe_sread_char failed (task)");

        if (task == COLOR_BOSS_YES) {
            rval = receive_prob (s, &probid, &ncount, &ecount, &elist);
            COLORcheck_rval (rval, "receive_prob failed");
        } else {
            printf ("No more work - shut down\n");
            goto CLEANUP;
        }

        COLORsafe_sclose (s);

        printf ("Received Graph %d: %d nodes, %d edges\n",
                 probid, ncount, ecount);
        fflush (stdout);

        printf ("Sleep 10 seconds, to simulate work\n"); fflush (stdout);
        sleep (10);

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

        printf ("Send boss the id (%d) of completed problem\n", probid);
        fflush (stdout);

        rval = COLORsafe_swrite_string (s, myname);
        COLORcheck_rval (rval, "COLORsafe_swrite_string failed (NAME)");

        rval = COLORsafe_swrite_char (s, COLOR_BOSS_RECEIVE);
        COLORcheck_rval (rval, "COLORsafe_swrite_char failed (RECEIVE)");

        rval = send_result (s, probid);
        COLORcheck_rval (rval, "send_result failed");

        COLORsafe_sclose (s);
    }


CLEANUP:

    COLOR_IFFREE (elist, int);
    return rval;
}

static int receive_prob (COLOR_SFILE *s, int *probid, int *ncount, int *ecount,
        int **pelist)
{
    int i, rval = 0;
    int *elist = (int *) NULL;

    rval = COLORsafe_sread_int (s, probid);
    COLORcheck_rval (rval, "COLORsafe_sread_int failed (probid)");

    rval = COLORsafe_sread_int (s, ncount);
    COLORcheck_rval (rval, "COLORsafe_sread_int failed (ncount)");

    rval = COLORsafe_sread_int (s, ecount);
    COLORcheck_rval (rval, "COLORsafe_sread_int failed (ecount)");

    elist = COLOR_SAFE_MALLOC (2*(*ecount), int);
    COLORcheck_NULL (elist, "out of memory for elist");

    for (i = 0; i < *ecount; i++) {
        rval = COLORsafe_sread_int (s, &(elist[2*i]));
        COLORcheck_rval (rval, "COLORsafe_sread_int failed (elist0)");

        rval = COLORsafe_sread_int (s, &(elist[2*i+1]));
        COLORcheck_rval (rval, "COLORsafe_sread_int failed (elist1)");
    }
    *pelist = elist;

CLEANUP:

    return rval;
}

static int send_result (COLOR_SFILE *s, int probid)
{
    static int rval = 0;

    rval = COLORsafe_swrite_int (s, probid);
    COLORcheck_rval (rval, "COLORsafe_swrite_int failed (probid)");

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
