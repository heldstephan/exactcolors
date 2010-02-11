#ifndef __COLOR_PRIVATE_H
#define __COLOR_PRIVATE_H

#include "bbsafe.h"
#include "lp.h"
#include "mwis.h"


typedef struct colordata colordata;
#define MAX_PNAME_LEN 256

struct colordata {

   /* The id of the node in the B&B tree.*/
   int id;
   int depth;

   enum {
      initialized        = 0,
      LP_bound_estimated = 1,
      LP_bound_computed  = 2,
      finished           = 3,
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

void init_colordata(colordata* cd);

int build_lp(colordata* cd);

void free_lbcolordata(colordata* cd);

void free_colordata(colordata* cd);

int create_branches(colordata* cd);

int send_colordata(COLOR_SFILE *s, colordata* cd, int include_best);

int receive_colordata(COLOR_SFILE *s, colordata* cd, 
                      int adopt_id,int include_best);

int compute_coloring(colordata* root_cd,int parallel);

int print_colors(COLORset* cclasses, int ccount);

void COLORset_dbg_lvl(int dbglvl);

void COLORset_cclasses_outfile(char* outfile);

void COLORset_write_mwis(int write_mwis);


#endif
