#include <stddef.h>
#include <limits.h>
#include <string.h>

#include "color_defs.h"

#include "color_parms.h"

void COLORparms_init(COLORparms* parms)
{
   parms->write_mwis                = 0;
   parms->branch_with_same_sequence = 0;
   parms->initial_upper_bound       = INT_MAX;
   parms->parallel_branching        = 0;


   parms->edgefile                  = (char*) NULL;    
   parms->outfile                   = (char*) NULL;    
   parms->cclasses_infile           = (char*) NULL;    
   parms->cclasses_outfile           = (char*) NULL;    
   parms->color_infile              = (char*) NULL;    
   parms->backupdir                 = (char*) NULL;       
}

void COLORparms_free(COLORparms* parms)
{
   COLOR_IFFREE(parms->edgefile,char);
   COLOR_IFFREE(parms->outfile,char);   
   COLOR_IFFREE(parms->cclasses_infile,char);
   COLOR_IFFREE(parms->cclasses_outfile,char);
   COLOR_IFFREE(parms->color_infile,char);
   COLOR_IFFREE(parms->backupdir,char);
}

static
int copy_string(char** dst, const char* src)
{
   int rval = 0;
   int len = strlen(src) + 1;
   COLOR_IFFREE(*dst,char);
   *dst = (char*) COLOR_SAFE_MALLOC(len,char);
   COLORcheck_NULL(*dst,"Failed to allocate dst");
   strcpy(*dst,src);
 CLEANUP:
   return rval;
}

int COLORparms_set_outfile(COLORparms* parms,const char* filename)
{
   return copy_string(&(parms->outfile),filename);
}

int COLORparms_set_edgefile(COLORparms* parms,const char* filename)
{
   return copy_string(&(parms->edgefile),filename);
}

int COLORparms_set_cclasses_infile(COLORparms* parms,const char* filename)
{
   return copy_string(&(parms->cclasses_infile),filename);
}

int COLORparms_set_cclasses_outfile(COLORparms* parms,const char* filename)
{
   return copy_string(&(parms->cclasses_outfile),filename);
}

int COLORparms_set_color_infile(COLORparms* parms,const char* filename)
{
   return copy_string(&(parms->color_infile),filename);
}

int COLORparms_set_backupdir(COLORparms* parms,const char* filename)
{
   return copy_string(&(parms->backupdir),filename);
}
int COLORparms_set_initial_upper_bound(COLORparms* parms,int bound)
{
   parms->initial_upper_bound = bound;
   return 0;
}

int COLORparms_set_write_mwis(COLORparms* parms,int write_mwis)
{
   parms->write_mwis = write_mwis;
   return 0;
}
