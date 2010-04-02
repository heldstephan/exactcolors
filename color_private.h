#ifndef __COLOR_PRIVATE_H
#define __COLOR_PRIVATE_H
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

#include "bbsafe.h"
#include "lp.h"
#include "mwis.h"

#include "heap.h"
#include "color_parms.h"

typedef struct colordata colordata;
#define MAX_PNAME_LEN 128

struct colordata {

   /* The id of the node in the B&B tree.*/
   int id;
   int depth;

   enum {
      initialized              = 0,
      LP_bound_estimated       = 1,
      LP_bound_computed        = 2,
      submitted_for_branching  = 3,
      finished                 = 4,
   } status;


   /* The instance graph */
   int ncount;
   int ecount;
   int *elist;
   int *orig_node_ids;


   /* The column generation LP. */
   COLORlp * lp;
   double *coef;
   double *pi;

   COLORNWT  mwis_pi_scalef;
   COLORNWT *mwis_pi;
   COLORNWT  lower_bound;
   COLORNWT  upper_bound;
   COLORNWT  lower_scaled_bound;
   double    dbl_safe_lower_bound;
   double    dbl_est_lower_bound;

   /* The MWIS instances. */
   MWISenv*  mwis_env;
   int       ccount;
   COLORset *cclasses;
   int       gallocated;
   COLORset *newsets;
   int       nnewsets;

   COLORset *bestcolors;
   int       nbestcolors;


   const COLORset *debugcolors;
   int             ndebugcolors;
   int   opt_track;

   int maxiterations;
   int retirementage;

   int v1,v2;
   colordata*       parent;
   colordata*       same_children;
   int              nsame;
   colordata*       diff_children;
   int              ndiff;

   char             pname[MAX_PNAME_LEN];
};

struct COLORproblem {
   COLORparms    parms;
   int           ncolordata;
   colordata     root_cd;
   COLORNWTHeap* br_heap;
   double        key_mult;
   COLORNWT      global_upper_bound;
};

void init_colordata(colordata* cd);

int set_id_and_name(colordata*  cd, 
                    int         id,
                    const char* name);

int init_unique_colordata(colordata*  cd, 
                          int         id,
                          const char* name);

int build_lp(colordata* cd);

void free_lbcolordata(colordata* cd);
void free_children_data(colordata* cd);
void free_colordata(colordata* cd);

int create_branches(colordata* cd, COLORproblem* problem);

int send_colordata(COLOR_SFILE *s, colordata* cd, int include_best);

int receive_colordata(COLOR_SFILE *s, colordata* cd,
                      int adopt_id,int include_best,
                      COLORproblem* problem);

int compute_coloring(COLORproblem* problem);

int print_colors(COLORset* cclasses, int ccount);

void COLORset_cclasses_outfile(char* outfile);

void COLORset_backupdir(char* backupdir);

const char* COLORget_backupdir(void );

int recover_colordata(colordata* cd,COLORproblem* problem);

int backup_colordata(colordata* cd, COLORproblem* problem);

int write_root_LP_snapshot(colordata* cd, COLORparms* parms, int add_timestamp);


void COLORset_write_mwis(int write_mwis);


#endif
