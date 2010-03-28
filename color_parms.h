#ifndef __COLOR_PARMS_H
#define __COLOR_PARMS_H
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

enum COLORBranchingStrategy {
   COLOR_min_strategy    = 0,
   COLOR_no_branching    = COLOR_min_strategy,

   COLOR_min_lb_strategy = 1,
   COLOR_dfs_strategy    = 2,
   COLOR_hybrid_strategy = 3,
   COLOR_max_strategy    = 4
};

typedef struct COLORparms {
   int      write_mwis;
   int      initial_upper_bound;
   int      parallel_branching;
   int      branch_with_same_sequence;
   int      branching_strategy;

   int      upper_bounds_only;

   double   branching_cpu_limit;

   char *edgefile;
   char *outfile;
   char *cclasses_infile;
   char *cclasses_outfile;
   char *color_infile;
   char *backupdir;



} COLORparms;

void COLORparms_init(COLORparms*);
void COLORparms_free(COLORparms*);

int COLORparms_set_outfile(COLORparms* parms,         const char* filename);
int COLORparms_set_edgefile(COLORparms* parms,        const char* filename);
int COLORparms_set_cclasses_infile(COLORparms* parms, const char* filename);
int COLORparms_set_cclasses_outfile(COLORparms* parms,const char* filename);
int COLORparms_set_color_infile(COLORparms* parms,    const char* filename);
int COLORparms_set_backupdir(COLORparms* parms,       const char* filename);


int COLORparms_set_initial_upper_bound(COLORparms* parms,int bound);
int COLORparms_set_write_mwis(COLORparms* parms,int write_mwis);
int COLORparms_set_parallel(COLORparms* parms,int parallel);
int COLORparms_set_branching_cpu_limit(COLORparms* parms, double branching_cpu_limit);
int COLORparms_set_branching_strategy(COLORparms* parms, int strategy);


#endif
