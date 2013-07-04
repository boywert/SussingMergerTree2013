#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


long total_halos(char filename[MAXSTRING]);
int lost_particles(long prog_halo, long merged_halo, int step);


long total_halos(char filename[MAXSTRING])
{
  FILE *fp;
  long totalhalo;
  char line[MAXSTRING];
  fp = fopen(filename,"r");
  if(fp == NULL)
    {
      fprintf(stderr,"could not open file %s\nexiting!\n",filename);
      exit(0);
    }
  fgets(line,MAXSTRING,fp);
  fclose(fp);
  sscanf(line, "%ld",  &totalhalo);

  return totalhalo;
 
}


int read_main_progenitor(char filename[MAXSTRING],long totalhalo)
{
  FILE *fpin;
  char line[MAXSTRING];
  long i,total_line;
  long item,item_parent;
  
  total_line = 0;
  
  for(i=0;i<totalhalo;i++)
    {
      main_progeniter[i] = -1;
    }
  fprintf(stderr,"reading file %s ...\n",filename);

  fpin = fopen(filename, "r");
  if(fpin == NULL)
    {
      fprintf(stderr,"could not open file %s\nexiting!\n",filename);
      exit(0);
    }
  while((fgets(line,MAXSTRING,fpin)) != NULL)
    {
      sscanf(line,"%ld %ld",&item, &item_parent);
      if(item < totalhalo)
	{
	  //printf("%ld %ld\n",item,totalhalo);
	  main_progeniter[item] = item_parent;
	}
      //main_progeniter[item] = item_parent;
      //printf("%ld %ld %ld\n",item,main_progeniter[item],totalhalo);
    }
  
  fclose(fpin);
  return 0;
}


int read_particles(char filename[MAXSTRING], int isimu)
{
  FILE *fpin;
  char  line[MAXSTRING];
  long  nPartInHalo, nPartInUse, ipart, jpart, ihalo, Pid, Ptype, numGoodHalos, haloid;
  long  PidMin_local=1234567890;
  long  PidMax_local=-1;
  
  fprintf(stderr,"o reading file %s ...",filename);
  
  fpin = fopen(filename,"r");
  if(fpin == NULL)
    {
      fprintf(stderr,"could not open file %s\nexiting!\n",filename);
      exit(0);
    }
  
  /* reset all variables */
  nHalos[isimu] = 0;
  ihalo         = -1;
  halos[isimu]   = NULL;
  
  
  /* get the first line from file */
  fgets(line,MAXSTRING,fpin);
  
  /* for AHF_particles files the first line is numGoodHalos which we can happily ignore */
  if(strncmp(line,"#",1) != 0 && sscanf(line,"%ld %ld",&haloid,&nPartInHalo) == 1)
    fgets(line,MAXSTRING,fpin);  
  
  do {
    if(strncmp(line,"#",1) != 0)
      {
      
	/* has a haloid been written before the number of particles (true for GalaxiesGoingMAD)*/
	if(sscanf(line,"%ld %ld",&haloid,&nPartInHalo) == 1)
	  {
	    /* if not, just get the number of particles */
	    sscanf(line,"%ld",&nPartInHalo);
	  }
      
	/* found yet another halo */
	ihalo++;
	nHalos[isimu] += 1;
	halos[isimu]   = (HALOptr) realloc(halos[isimu], nHalos[isimu]*sizeof(HALOS));
      
	/* store npart and allocate Pid[] accordingly */
	halos[isimu][ihalo].Pid   = (long *) calloc(nPartInHalo, sizeof(long));
      
	/* read all their id's */
	nPartInUse = 0;
	for(ipart=0; ipart<nPartInHalo; ipart++)
	  {
	    /* read line containing ID and possibly some more information */
	    fgets(line,MAXSTRING,fpin);
        
	    /* check whether we are able to read the particle type, too */
	    if(sscanf(line,"%ld %ld",&Pid,&Ptype) == 1)
	      {
		/* if not, set Ptype to some random number to be ignored */
		sscanf(line,"%ld",&Pid);
		Ptype = -1;
	      }
        
	    // here we can restrict the cross-correlation to a ceratain sub-set of all particles
	    if(Ptype == 1)
	      {  
		halos[isimu][ihalo].Pid             = (long *) realloc(halos[isimu][ihalo].Pid, (nPartInUse+1)*sizeof(long));
		halos[isimu][ihalo].Pid[nPartInUse] = Pid;
           
		if(abs(halos[isimu][ihalo].Pid[nPartInUse]) > PidMax) PidMax = abs(halos[isimu][ihalo].Pid[nPartInUse]);
		if(abs(halos[isimu][ihalo].Pid[nPartInUse]) < PidMin) PidMin = abs(halos[isimu][ihalo].Pid[nPartInUse]);
		if(abs(halos[isimu][ihalo].Pid[nPartInUse]) > PidMax_local) PidMax_local = abs(halos[isimu][ihalo].Pid[nPartInUse]);
		if(abs(halos[isimu][ihalo].Pid[nPartInUse]) < PidMin_local) PidMin_local = abs(halos[isimu][ihalo].Pid[nPartInUse]);
           
		nPartInUse++;
	      }
	  }
      
	/* store number of particles in halo */
	halos[isimu][ihalo].npart = nPartInUse;
      
#ifdef DEBUG
	fprintf(stderr,"  => halo %ld of simu #%ld contains %ld dm+gas particles (expected %ld particles in total)\n",ihalo,isimu,nPartInUse,nPartInHalo);
#endif
      }
  } while( fgets(line,MAXSTRING,fpin) != NULL);
  
  fclose(fpin);
  
  fprintf(stderr," done (full ID range = %ld -> %ld, local ID range = %ld -> %ld)\n",PidMin,PidMax,PidMin_local,PidMax_local);
  
  return(1);
}



/*==================================================================================================
 * particle_halo_mapping:
 *
 *  for each particle remember to which halo(s) it belongs
 *
 *==================================================================================================*/
int particle_halo_mapping(int isimu)
{
  long ihalo, ipart, jpart;
  
  fprintf(stderr,"o creating particle<->halo mapping for file %d ...",isimu);
  
  
  /* recording every halo it belongs to */
  for(ihalo=0; ihalo<nHalos[isimu]; ihalo++)
    {
      for(jpart=0; jpart<halos[isimu][ihalo].npart; jpart++)
	{
	  ipart = abs(halos[isimu][ihalo].Pid[jpart]);
      
	  parts[isimu][ipart].nhalos++;
	  parts[isimu][ipart].Hid = (long *) realloc(parts[isimu][ipart].Hid, parts[isimu][ipart].nhalos*sizeof(long));
      
	  parts[isimu][ipart].Hid[parts[isimu][ipart].nhalos-1] = ihalo;
	}
    }
  
  fprintf(stderr," done\n");
  return(1);
}

/*==================================================================================================
 * cross_correlation:
 *
 *  for each halo at isimu=0 figure out how many particles are in common with khalo at isimu=1
 *
 *==================================================================================================*/

int lost_particles(long prog_halo, long merged_halo, int step)
{
  long  khalo, ipart, jpart;
  long *common;
  FILE *fpout;
  char outname[MAXSTRING];
  long commonpart;
  long i;
  
  /* common records how many particles ihalo(isimu=0) has in common with khalo(isimu=1) */
  commonpart = 0;
  for(jpart=0; jpart<halos[1][merged_halo].npart; jpart++)
    {
      //printf("%ld < %ld \n",jpart, halos[1][ihalo].npart);
      ipart = halos[1][merged_halo].Pid[jpart];
      for(i=0;i < parts[0][ipart].nhalos; i++)
	{
	  if(parts[0][ipart].Hid[i] == prog_halo)
	    {
	      commonpart++;
	    }
	}
      /* ipart belongs to nhalos halos in isimu=1 */
	             
      //khalo          = parts[0][ipart].Hid[jhalo];
      //if(jhalo == khalo)
      //	commonpart += 1;
    }
  if(commonpart > 0)
    {
      //printf("%ld => %ld : %ld from %ld => %f\n",prog_halo,merged_halo,commonpart,halos[0][prog_halo].npart, (double) (halos[0][prog_halo].npart - commonpart) / (double) halos[0][prog_halo].npart );
      delta_out[step] += (double) (halos[0][prog_halo].npart - commonpart) / (double) halos[0][prog_halo].npart;
      merging_rate[step] +=  (double) (halos[1][merged_halo].npart - halos[0][prog_halo].npart) / (double) halos[0][prog_halo].npart;
      return 1;
    }
  else
    {
      return 0;
    }
  /* write info to file 
  for(khalo=0; khalo<nHalos[1]; khalo++)
    {
      if(common[khalo] > MINCOMMON)
	fprintf(fpout,"%12ld  %12ld  %12ld  %12ld  %12ld\n",
		ihalo,
		halos[0][ihalo].npart,
		common[khalo],
		khalo,
		halos[1][khalo].npart);
    }
    */
}

