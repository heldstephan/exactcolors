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

static char *edgefile = (char *) NULL;
static int debug = 0;

int main (int ac, char **av);
static int send_prob (COLOR_SFILE *s, int current, int ncount, int ecount,
    int *elist);
static int receive_result (COLOR_SFILE *s, int *probid);
static int parseargs (int ac, char **av);
static void usage (char *f);

int main (int ac, char **av)
{
    int rval = 0;
    int ncount, ecount;
    int *elist = (int *) NULL;
    COLOR_SPORT *lport = (COLOR_SPORT *) NULL;
    COLOR_SFILE *s;
    char request, grunt[1024];
    int donecount = 0, docount = 5;
    int probid = 0, current = 0;

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    rval = COLORread_dimacs (edgefile, &ncount, &ecount, &elist,
                             (int **) NULL);
    COLORcheck_rval (rval, "COLORread_dimacs failed");

    printf ("\nTEST BOSS CODE\n\n");
    fflush (stdout);

    lport = COLORsafe_snet_listen (COLOR_BOSS_PORT);
    if (lport == (COLOR_SPORT *) NULL) {
        fprintf (stderr, "COLORsafe_snet_listen failed\n");
        rval = 1; goto CLEANUP;
    }

    printf ("Number of Test Problems: %d\n", docount); fflush (stdout);

    do {
        s = COLORsafe_snet_receive (lport);
        if (!s) {
            fprintf (stderr, "COLORsafe_snet_receive failed, ignoring\n");
            continue;
        }

        if (COLORsafe_sread_string (s, grunt, 1023)) {
            fprintf (stderr, "COLORsafe_sread_char string, abort con\n");
            COLORsafe_sclose (s);
            continue;
        }

        if (COLORsafe_sread_char (s, &request)) {
            fprintf (stderr, "COLORsafe_sread_char failed, abort con\n");
            COLORsafe_sclose (s);
            continue;
        }

        switch (request) {
        case COLOR_BOSS_SEND:
            if (current < docount) {
                rval = send_prob (s, current, ncount, ecount, elist);
                if (rval) {
                    fprintf (stderr, "send_cut failed - abort con\n");
                } else {
                    printf ("Sent test %d to %s (remaining %d)\n",
                         current, grunt, docount - donecount);
                    fflush (stdout);
                    current++;
                }
            } else {
                rval = COLORsafe_swrite_char (s, COLOR_BOSS_NO);
                if (rval) {
                    fprintf (stderr, "BOSS_NONE write failed - abort con\n");
                }
            }
            COLORsafe_sclose (s);
            break;
        case COLOR_BOSS_RECEIVE:
            rval = receive_result (s, &probid);
            if (rval) {
                fprintf (stderr, "receive_result failed - abort connection\n");
                COLORsafe_sclose (s);
                break;
            } else {
                printf ("Received result %d from %s\n", probid, grunt);
                fflush (stdout);
                COLORsafe_sclose (s);
                donecount++;
            }
            break;
        case COLOR_BOSS_EXIT:
            printf ("Shutting down the test boss\n"); fflush (stdout);
            COLORsafe_sclose (s);
            goto CLEANUP;
        default:
            fprintf (stderr, "Invalid request %c\n", request);
        }
    } while (donecount != docount);

CLEANUP:

    COLORsafe_snet_unlisten (lport);
    COLOR_IFFREE (elist, int);
    return rval;
}

static int send_prob (COLOR_SFILE *s, int current, int ncount, int ecount,
        int *elist)
{
    int i, rval = 0;

    rval = COLORsafe_swrite_char (s, COLOR_BOSS_YES);
    COLORcheck_rval (rval, "COLORsafe_swrite_int failed (YES)");

    rval = COLORsafe_swrite_int (s, current);
    COLORcheck_rval (rval, "COLORsafe_swrite_int failed (cut id)");

    rval = COLORsafe_swrite_int (s, ncount);
    COLORcheck_rval (rval, "COLORsafe_swrite_int failed (ncount)");

    rval = COLORsafe_swrite_int (s, ecount);
    COLORcheck_rval (rval, "COLORsafe_swrite_int failed (ecount)");

    for (i = 0; i < ecount; i++) {
        rval = COLORsafe_swrite_int (s, elist[2*i]);
        COLORcheck_rval (rval, "COLORsafe_swrite_int failed (ecount0)");

        rval = COLORsafe_swrite_int (s, elist[2*i+1]);
        COLORcheck_rval (rval, "COLORsafe_swrite_int failed (ecount1)");
    }

CLEANUP:

    return rval;
}

static int receive_result (COLOR_SFILE *s, int *probid)
{
    int rval = 0;

    rval = COLORsafe_sread_int (s, probid);
    COLORcheck_rval (rval, "COLORsafe_sread_int failed (probid)");

CLEANUP:

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int rval = 0;

    while ((c = getopt (ac, av, "d")) != EOF) {
        switch (c) {
        case 'd':
           /* each -d increases the verbosity by one.*/
           ++debug;
            break;
        default:
            rval = 1;
            goto CLEANUP;
        }
    }

    if (ac <= optind) {
        rval = 1; goto CLEANUP;
    } else {
        edgefile = av[optind++];
    }

CLEANUP:

    if (rval) { usage (av[0]); }
    return rval;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage %s: [-see below-] edge_file\n", f);
    fprintf (stderr, "   -d    turn on debugging\n");
}

int COLORdbg_lvl() {
   return debug;
}
