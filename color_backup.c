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

#include <sys/stat.h>
#include <time.h>

#include "color_private.h"

static char* backupdir = 0;


void COLORset_backupdir(char* dir)
{
   backupdir = dir;
}

const char* COLORget_backupdir(void )
{
   return backupdir;
}


COLOR_MAYBE_UNUSED
static int write_colordata_to_file(colordata* cd, 
				   FILE* file) {
   int prval = 0;
   int rval  = 0;
   int i;
   

   prval = fprintf(file,"id %d\n",cd->id);
   COLORcheck_fileio(prval,"Failed in fprintf");

   if (cd->parent) {
      prval = fprintf(file,"parent_id %d\n",cd->parent->id);
      COLORcheck_fileio(prval,"Failed in fprintf");

   } else {
      prval = fprintf(file,"parent_id -1\n");
      COLORcheck_fileio(prval,"Failed in fprintf");

   }
   
   prval = fprintf(file,"pname %s\n",cd->pname);
   COLORcheck_fileio(prval,"Failed in fprintf");


   prval = fprintf(file,"depth %d\n",cd->depth);
   COLORcheck_fileio(prval,"Failed in fprintf");


   prval = fprintf(file,"status %u\n",cd->status);
   COLORcheck_fileio(prval,"Failed in fprintf");

   prval = fprintf(file,"ecount %d\n",cd->ecount);
   COLORcheck_fileio(prval,"Failed in fprintf");


   for (i = 0; i < cd->ecount; ++ i) {
      prval = fprintf(file,"e %d %d\n",cd->elist[2*i], cd->elist[2*i+1]);
      COLORcheck_fileio(prval,"Failed in fprintf");
   }

   prval = fprintf(file,"ncount %d\n",cd->ncount);
   COLORcheck_fileio(prval,"Failed in fprintf");

   prval = fprintf(file,"orig_node_ids ");
   COLORcheck_fileio(prval,"Failed in fprintf");
   for (i = 0; i < cd->ncount; ++ i) {
      prval = fprintf(file," %d",cd->orig_node_ids[i]);
      COLORcheck_fileio(prval,"Failed in fprintf");
   }
   prval = fprintf(file,"\n");
   COLORcheck_fileio(prval,"Failed in fprintf");


   prval = fprintf(file,"lower_bound %d\n",cd->lower_bound);
   COLORcheck_fileio(prval,"Failed in fprintf");


   prval = fprintf(file,"upper_bound %d\n",cd->upper_bound);
   COLORcheck_fileio(prval,"Failed in fprintf");


   prval = fprintf(file,"dbl_safe_lower_bound %f\n",cd->dbl_safe_lower_bound);
   COLORcheck_fileio(prval,"Failed in fprintf");

   prval = fprintf(file,"dbl_est_lower_bound %f\n",cd->dbl_est_lower_bound);
   COLORcheck_fileio(prval,"Failed in fprintf");

   prval = fprintf(file,"ccount %d\n",cd->ccount);
   COLORcheck_fileio(prval,"Failed in fprintf");

   for (i = 0; i < cd->ccount; ++ i) {
      int j;
      prval = fprintf(file,"count %d\n",cd->cclasses[i].count);
      COLORcheck_fileio(prval,"Failed in fprintf");

      prval = fprintf(file,"members");
      COLORcheck_fileio(prval,"Failed in fprintf");

      for (j = 0; j < cd->cclasses[i].count; ++j) {
	 prval = fprintf(file," %d",cd->cclasses[i].members[j]);
	 COLORcheck_fileio(prval,"Failed in fprintf");
      }
      prval = fprintf(file,"\n");
      COLORcheck_fileio(prval,"Failed in fprintf");
   }
   
   
   prval = fprintf(file,"nbestcolors %d\n",cd->nbestcolors);
   COLORcheck_fileio(prval,"Failed in fprintf");
   
   for (i = 0; i < cd->nbestcolors; ++ i) {
      int j;
      prval = fprintf(file,"bestcolors_count %d\n",cd->bestcolors[i].count);
      COLORcheck_fileio(prval,"Failed in fprintf");
      prval = fprintf(file,"bestcolors_members");
      COLORcheck_fileio(prval,"Failed in fprintf");
      for (j = 0; j < cd->bestcolors[i].count; ++j) {
	 prval = fprintf(file," %d",cd->bestcolors[i].members[j]);
	 COLORcheck_fileio(prval,"Failed in fprintf");
      }
      prval = fprintf(file,"\n");
      COLORcheck_fileio(prval,"Failed in fprintf");
   }
   
   prval = fprintf(file,"v1 %d\n",cd->v1);
   COLORcheck_fileio(prval,"Failed in fprintf");
   prval = fprintf(file,"v2 %d\n",cd->v2);
   COLORcheck_fileio(prval,"Failed in fprintf");

   prval = fprintf(file,"nsame %d\n",cd->nsame);
   COLORcheck_fileio(prval,"Failed in fprintf");
   
   prval = fprintf(file,"same_children");
   COLORcheck_fileio(prval,"Failed in fprintf");
   for (i = 0; i < cd->nsame; ++ i) {
      prval = fprintf(file," %d",cd->same_children[i].id);
      COLORcheck_fileio(prval,"Failed in fprintf");
   }
   prval = fprintf(file,"\n");
   COLORcheck_fileio(prval,"Failed in fprintf");

   prval = fprintf(file,"ndiff %d\n",cd->ndiff);
   COLORcheck_fileio(prval,"Failed in fprintf");

   prval = fprintf(file,"diff_children");
   COLORcheck_fileio(prval,"Failed in fprintf");
   for (i = 0; i < cd->ndiff; ++ i) {
      prval = fprintf(file," %d\n",cd->diff_children[i].id);
      COLORcheck_fileio(prval,"Failed in fprintf");
   }
   /* prval = fprintf(file,"\n"); */
   /* COLORcheck_fileio(prval,"Failed in fprintf"); */

 CLEANUP:
   return rval;
}

int backup_colordata(colordata* cd) {
   
   char bkp_filename[256] = "";

   int   prval = 0;
   int   rval  = 0;
   FILE* file  = (FILE*) NULL;
   
   if (backupdir) {
      char filename[256];
      
      sprintf(bkp_filename,"%s",cd->pname);

      if(!COLORdir_exists(backupdir)) {
	 prval = mkdir(backupdir,(S_IRUSR | S_IWUSR | S_IXUSR));
	 COLORcheck_fileio(prval,"Failed to mkdir");
      }
    
      /* prval = sprintf(filename,"gzip >|%s/%s.%d.gz", */
      /* 		      backupdir, cd->pname, cd->id); */
      /* COLORcheck_fileio(prval,"Failed in sprintf");	  */
      /* file = popen(filename,"w"); */
      /* COLORcheck_NULL(file,"Failed to open file"); */

      prval = sprintf(filename,"%s/%s.%d",
      		      backupdir, cd->pname, cd->id);
      COLORcheck_fileio(prval,"Failed in sprintf");
      if(COLORfile_exists(filename)) {
         prval = sprintf(bkp_filename,"%s.bkp",
                         filename);
         COLORcheck_fileio(prval,"Failed in sprintf");
      
         printf("Renaming %s to %s before update.\n",
                filename, bkp_filename);
         rval = rename(filename,bkp_filename);
         COLORcheck_rval(rval, "Failed to rename file.");
      }
      file = fopen(filename,"w"); 

      rval = write_colordata_to_file
	 (cd,file);
      COLORcheck_rval(rval,"Failed in write_colordata");
   }
CLEANUP:
   if (file) {
      /* pclose(file); */
      fclose(file);
   }

   if(!rval && strcmp(bkp_filename,cd->pname) && 
      COLORfile_exists(bkp_filename)) {
      remove(bkp_filename);
   }

   return rval;
}


COLOR_MAYBE_UNUSED
static int read_colordata_from_file(colordata* cd, 
                                    COLORproblem* problem,
                                    FILE* file) 
{
   const int LINE_SIZE = 32769;
   int rval = 0;
   int prval = 0;
   char buf[LINE_SIZE], *p;
   int  ecount_test = 0;
   int  cclass_i = 0;
   int  bestc_i = 0;
   const char* delim = " \t\n";

   while (fgets (buf,LINE_SIZE, file) != (char *) NULL) {
      char* data = (char *) NULL;

      p = buf;
      data = strtok(p,delim); /* get 'p' */

      if (!strcmp(data,"id")) {
         data = strtok((char*)NULL,delim);
         cd->id = atoi(data);
      }

      else if (!strcmp(data,"pname")) {
         data = strtok((char*)NULL,delim);
         prval = sprintf(cd->pname,data);
         COLORcheck_fileio(prval,"Failed in sprintf");
      }
      else if (!strcmp(data,"depth")) {
         data = strtok((char*)NULL,delim);
         cd->depth = atoi(data);
      }
      else if (!strcmp(data,"status")) {
         data = strtok((char*)NULL,delim);
         cd->status = atoi(data);
      }
      else if (!strcmp(data,"ecount")) {
         data = strtok((char*)NULL,delim);
         cd->ecount = atoi(data);
         COLOR_IFFREE(cd->elist,int);
         cd->elist = (int*) COLOR_SAFE_MALLOC (2 * cd->ecount, int);
      }

      else if (!strcmp(data,"e")) {
         data = strtok((char*)NULL,delim);
         cd->elist[2*ecount_test] = atoi(data);
         data = strtok((char*)NULL,delim);
         cd->elist[2*ecount_test + 1] = atoi(data);
         ecount_test++;
      }

      else if (!strcmp(data,"ncount")) {
         data = strtok((char*)NULL,delim);
         cd->ncount = atoi(data);
         COLOR_IFFREE(cd->orig_node_ids,int);
         cd->orig_node_ids = (int*) COLOR_SAFE_MALLOC (cd->ncount, int);
      }

      else if (!strcmp(data,"orig_node_ids")) {
         int i;
         for (i = 0; i < cd->ncount; ++i) {
            data = strtok((char*)NULL,delim);
            cd->orig_node_ids[i] = atoi(data);
         }
      }
      else if (!strcmp(data,"lower_bound")) {
         data = strtok((char*)NULL,delim);
         cd->lower_bound = atoi(data);
      }

      else if (!strcmp(data,"upper_bound")) {
         data = strtok((char*)NULL,delim);
         cd->upper_bound = atoi(data);
      }
      else if (!strcmp(data,"dbl_safe_lower_bound")) {
         data = strtok((char*)NULL,delim);
         cd->dbl_safe_lower_bound = atof(data);
      }
      else if (!strcmp(data,"dbl_est_lower_bound")) {
         data = strtok((char*)NULL,delim);
         cd->dbl_est_lower_bound = atof(data);
      }

      else if (!strcmp(data,"ccount")) {
         data = strtok((char*)NULL,delim);
         cd->ccount = atoi(data);
         COLOR_IFFREE(cd->cclasses,COLORset);
         if (cd->ccount) cd->cclasses = COLOR_SAFE_MALLOC(cd->ccount, COLORset);
      }
            
      else if (!strcmp(data,"count")) {
         data = strtok((char*)NULL,delim);
         cd->cclasses[cclass_i].count = atoi(data);
         if (cd->cclasses[cclass_i].count) cd->cclasses[cclass_i].members = COLOR_SAFE_MALLOC(cd->cclasses[cclass_i].count, int);
      }

      else if (!strcmp(data,"members")) {
         int i;
         for (i = 0; i < cd->cclasses[cclass_i].count; ++i) {
            data = strtok((char*)NULL,delim);
            cd->cclasses[cclass_i].members[i] = atoi(data);
         }
         cclass_i++;
      }

      else if (!strcmp(data,"nbestcolors")) {
         data = strtok((char*)NULL,delim);
         cd->nbestcolors = atoi(data);
         COLOR_IFFREE(cd->bestcolors, COLORset);
         if (cd->nbestcolors) cd->bestcolors = COLOR_SAFE_MALLOC(cd->nbestcolors, COLORset);
      }
            
      else if (!strcmp(data,"bestcolors_count")) {
         data = strtok((char*)NULL,delim);
         cd->bestcolors[bestc_i].count = atoi(data);
         if (cd->bestcolors[bestc_i].count) cd->bestcolors[bestc_i].members = COLOR_SAFE_MALLOC(cd->bestcolors[bestc_i].count, int);
      }

      else if (!strcmp(data,"bestcolors_members")) {
         int i;
         for (i = 0; i < cd->bestcolors[bestc_i].count; ++i) {
            data = strtok((char*)NULL,delim);
            cd->bestcolors[bestc_i].members[i] = atoi(data);
         }
         bestc_i++;
      }

      else if (!strcmp(data,"v1")) {
         data = strtok((char*)NULL,delim);
         cd->v1 = atoi(data);
      }

      else if (!strcmp(data,"v2")) {
         data = strtok((char*)NULL,delim);
         cd->v2 = atoi(data);
      }

      else if (!strcmp(data,"nsame")) {
         data = strtok((char*)NULL,delim);
         cd->nsame = atoi(data);
         COLOR_IFFREE(cd->same_children,colordata);
         if (cd->nsame) {
            int i;
            cd->same_children = COLOR_SAFE_MALLOC(cd->nsame,colordata);
            for (i = 0; i < cd->nsame; ++i) {
               init_colordata(cd->same_children + i);
            }
         }
      }

      else if (!strcmp(data,"same_children")) {
         int i;
         for (i = 0; i < cd->nsame; ++i) {
            data = strtok((char*)NULL,delim);
            cd->same_children[i].id = atoi(data);
            strcpy(cd->same_children[i].pname,cd->pname);
            
            rval = recover_colordata(cd->same_children + i,problem);
            if (rval) {
               free_children_data(cd); rval = 0;
               goto CLEANUP;
            }
            cd->same_children[i].parent = cd;
         }
      }

      else if (!strcmp(data,"ndiff")) {
         data = strtok((char*)NULL,delim);
         cd->ndiff = atoi(data);
         COLOR_IFFREE(cd->diff_children,colordata);
         if (cd->ndiff) {
            int i;
            cd->diff_children = COLOR_SAFE_MALLOC(cd->ndiff,colordata);
            for (i = 0; i < cd->ndiff; ++i) {
               init_colordata(cd->diff_children + i);
            }

         }
      }

      else if (!strcmp(data,"diff_children")) {
         int i;
         for (i = 0; i < cd->ndiff; ++i) {
            data = strtok((char*)NULL,delim);
            cd->diff_children[i].id = atoi(data);
            strcpy(cd->diff_children[i].pname,cd->pname);

            rval = recover_colordata(cd->diff_children + i,problem);
            if (rval) {
               free_children_data(cd); rval = 0;
               goto CLEANUP;
            }
            cd->diff_children[i].parent = cd;
         }
      }
   }
 CLEANUP:
   return rval;
}

int recover_colordata(colordata* cd,COLORproblem* problem) {
   int rval  = 0;
   int prval = 0;
   FILE* file = (FILE*) NULL;
   if (backupdir) {
      char filename[256];
      prval = sprintf(filename,"%s/%s.%d",
      		      backupdir, cd->pname, cd->id);
      COLORcheck_fileio(prval,"Failed in sprintf");
      file = fopen(filename,"r"); 
      COLORcheck_NULL(file, "Failed to fopen");

      init_colordata(cd);

      rval = read_colordata_from_file(cd, problem, file);
      COLORcheck_rval(rval,"Failed in read_colordata_from_file");
   }
 CLEANUP:
   if (file) {
      fclose(file);
   }
   return rval;
}


int write_root_LP_snapshot(colordata* cd, COLORparms* parms, int add_timestamp)
{
   int rval = 0;
   if (parms->cclasses_outfile != (char*) NULL) {
      char   fname[256];
      const colordata* root_cd = cd;

      while (root_cd->parent) {
         root_cd = root_cd->parent;
      }

      if (add_timestamp) {
         char   timestr[256];
         time_t t;
         t = time(NULL);

         strftime(timestr, 256, "%Y%m%d%H%M", localtime(&t));
         sprintf(fname,"%s.%s",parms->cclasses_outfile,timestr);
      } else {
         sprintf(fname,"%s",parms->cclasses_outfile);
      }
      rval = COLORstable_write_stable_sets(cd->cclasses,cd->ccount,cd->ncount,
                                           fname,root_cd->pname);
      COLORcheck_rval(rval,"Failed in COLORstable_write_stable_sets");
   }
 CLEANUP:
   return rval;
}
