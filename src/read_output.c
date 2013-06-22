#include "allvars.h"

struct outputRAW output;

MyIDtype TotNavatars;
MyIDtype *Avatar;

int readRawOutputs(char file[MAXSTRING])
{
  FILE *fp;
  MyIDtype i,j,line,dummynumber;
  char buffer[MAXSTRING];
  char formatstr[MAXSTRING];
  time_t clock;
  struct stat attrib;
  line = 0;
  strncpy(formatstr,"",sizeof(formatstr));
  fp = fopen(file, "r");
  if(fp == NULL)
    {
      fprintf(stderr,"Error: Cannot read output file\nExiting (%llu)...\n",line);
      //exit(0);
      return 1;
    }
  strcpy(output.filename, file);
  stat(file, &attrib);
  output.timestamp = gmtime(&(attrib.st_mtime));
  fprintf(stderr,"Filename : %s\n",output.filename);
  clock = mktime(output.timestamp);
  
  fprintf(stderr,"Last modified (UCT) : %s\n",ctime(&clock));
  line++;
  strncpy(buffer,"",sizeof(buffer));
  if((fgets(buffer,MAXSTRING,fp)) == NULL)
    {
      fprintf(stderr,"Error: Cannot read output file\nExiting (%llu)...\n",line);
      //exit(0);
      fclose(fp);
      return 1;
    }
  strcat(formatstr,"%f");
  if(sscanf(buffer,formatstr,&(output.outputFormat)) < 0.99)
    {
      fprintf(stderr,"Error: invalid data format\nExiting (%llu)...\n",line);
      //exit(0);
      fclose(fp);
      return 1;
    }
  if((int) output.outputFormat == 1)
    {
      //printf("formatID = %d\n",output.outputFormat);
      line++;
      strncpy(buffer,"",sizeof(buffer));
      fgets(output.Descriptions,MAXSTRING,fp);
      fprintf(stderr,"Description: %s\n",output.Descriptions);
      line++;
      strncpy(buffer,"",sizeof(buffer));
      if((fgets(buffer,MAXSTRING,fp)) == NULL)
	{
	  fprintf(stderr,"Error: Cannot read nHalos\nExiting (%llu)...\n",line);
	  //exit(0);
	  fclose(fp);
	  return 1;
	}
      strncpy(formatstr,"",sizeof(formatstr));
      strcat(formatstr,"%llu");
      if(sscanf(buffer,formatstr,&(output.nHalos)) < 1)
	{
	  fprintf(stderr,"Error: Invalid format\nExiting (%llu)...\n",line);
	  //exit(0);
	  fclose(fp);
	  return 1;
	}
      if(output.nHalos < 1)
	{
	  fprintf(stderr,"Error: Invalid nHalos = %llu\nExiting (%llu)...\n",output.nHalos,line);
	  //exit(0);
	  fclose(fp);
	  return 1;
	}
      /*
      if(output.FirstSnapshot < FIRSTSNAP)
	{
	  fprintf(stderr,"Error: Invalid number of first snapshot = %llu\nExiting (%llu)...\n",output.FirstSnapshot,line);
	  //exit(0);
	  fclose(fp);
	  return 1;
	}
      if(output.LastSnapshot > LASTSNAP)
	{
	  fprintf(stderr,"Error: Invalid number of last snapshot = %llu\nExiting (%llu)...\n",output.LastSnapshot,line);
	  //exit(0);
	  fclose(fp);
	  return 1;
	}
      if(output.TotNsnapshots > NSNAPS || output.TotNsnapshots < 1 )
	{
	  fprintf(stderr,"Error: Invalid number of snapshot = %llu\nExiting (%llu)...\n",output.TotNsnapshots,line);
	  //exit(0);
	  fclose(fp);
	  return 1;
	}
      */
      //printf("nHalo : %llu\n",output.nHalos);
      (void) makeIDmap();

      if(output.nHalos != TotNhalos )
	{
	  fprintf(stderr,"Error: Invalid output nHalos= %llu != %llu\nExiting (%llu)...\n",output.nHalos,TotNhalosUsed,line);
	  //exit(0);
	  fclose(fp);
	  return 1;
	}
      
      output.progs = calloc(output.nHalos, sizeof(struct progtable));
      for(i=0;i<output.nHalos;i++)
	{
	  line++;
	  strncpy(buffer,"",sizeof(buffer));
	  if((fgets(buffer,MAXSTRING,fp)) == NULL)
	    {
	      fprintf(stderr,"Error: nHalos mismatched\nExiting (%llu)...\n",line);
	      //exit(0);
	      fclose(fp);
	      return 1;
	    }
	  strncpy(formatstr,"",sizeof(formatstr));
	  strcat(formatstr,"%llu %llu");
	  if ( (sscanf(buffer,formatstr, &(output.progs[i].haloID), &(output.progs[i].nProgs))) < 1)
	    {
	      fprintf(stderr,"Error: invalid format\nExiting (%llu)...\n",line);
	      //exit(0);
	      fclose(fp);
	      return 1;
	    }
	  dummynumber = IDsearch(output.progs[i].haloID);
	  if(dummynumber==NULLPOINT )
	    {
	      fprintf(stderr,"Error: invalid HaloID %llu\nExiting (%llu)...\n",output.progs[i].haloID,line);
	      //exit(0);
	      fclose(fp);
	      return 1;
	    }
	  else
	    {
	      output.progs[i].haloID = dummynumber;
	    }
	  //printf("HaloID: %llu has %llu progenitors\n",output.progs[i].haloID,output.progs[i].nProgs);
	  if(output.progs[i].nProgs > 0)
	    {
	      output.progs[i].progID = calloc(output.progs[i].nProgs, sizeof(MyIDtype));
	      for(j=0;j<output.progs[i].nProgs;j++)
		{
		  line++;
		  strncpy(buffer,"",sizeof(buffer));
		  if((fgets(buffer,MAXSTRING,fp)) == NULL)
		    {
		      fprintf(stderr,"Error: nProgs mismatched\nExiting (%llu)...\n",line);
		      //exit(0);
		      fclose(fp);
		      return 1;
		    }
		  strncpy(formatstr,"",sizeof(formatstr));
		  strcat(formatstr,"%llu");
		  if(sscanf(buffer,formatstr,&(output.progs[i].progID[j])) < 1)
		    {
		      fprintf(stderr,"Error: invalid format\nExiting (%llu)...\n",line);
		      //exit(0);
		      fclose(fp);
		      return 1;
		    } 
		  dummynumber = IDsearch(output.progs[i].progID[j]);
		  if(dummynumber == NULLPOINT)
		    {
		      fprintf(stderr,"Error: invalid HaloID %llu\nExiting (%llu)...\n",output.progs[i].progID[j],line);
		      //exit(0);
		      fclose(fp);
		      return 1;
		    }
		  else
		    {
		      output.progs[i].progID[j] = dummynumber;
		    }
		}
	    }

	}
      line++;
      strncpy(buffer,"",sizeof(buffer));
      fgets(buffer,MAXSTRING,fp);
      if(strstr(buffer, "END") == NULL)
	{
	  fprintf(stderr,"Error: line mismatched\nExiting (%llu)...\n",line);
	  //exit(0);
	  fclose(fp);
	  return 1;
	}
    }
  else
    {
      fprintf(stderr,"Error: unsupport formatID = %f\nExiting (%llu)...\n",output.outputFormat,line);
      //exit(0);
      fclose(fp);
      return 1;
    }
  fclose(fp); 
  return 0;
}
MyIDtype assignAvatarMain()
{
  MyIDtype i,j,loop1,loop2,avatardup;
  avatardup = 0;
  for(i=0;i<output.nHalos;i++)
    {
      for(j=0;j<output.progs[i].nProgs;j++)
	{
	  if(HaloTable[output.progs[i].progID[j]].ProgAvatarFlag == 0)
	    {
	      HaloTable[output.progs[i].progID[j]].ProgAvatarFlag = 1;
	    }
	  else if(HaloTable[output.progs[i].progID[j]].ProgAvatarFlag == 1)
	    {
	      HaloTable[output.progs[i].progID[j]].nAvatars++;
	      TotNavatars++;
	      HaloTable[output.progs[i].progID[j]].AvatarList = realloc(HaloTable[output.progs[i].progID[j]].AvatarList, HaloTable[output.progs[i].progID[j]].nAvatars*sizeof(MyIDtype));
	      Avatar = realloc(Avatar, TotNavatars*sizeof(MyIDtype));
	      HaloTable[output.progs[i].progID[j]].AvatarList[HaloTable[output.progs[i].progID[j]].nAvatars-1]=TotNavatars-1;
	      Avatar[TotNavatars-1] = HaloTable[output.progs[i].progID[j]].AvatarList[0];
	    }
	  else
	    {
	      exit(0);
	    }
	  break;
	}
    }
  for(i=0;i<TotNhalos;i++)
    {
      if(HaloTable[i].nAvatars > 1)
      	{
	  avatardup++;
      // printf("%f\t%llu\n",HaloTable[i].Mvir,HaloTable[i].nAvatars);
	}
    }
  printf("Main: TotNhalos %llu, TotNavatar = %llu,  ExtraAvatars = %llu, MultipleAvatars = %llu\n",TotNhalos,TotNavatars,TotNavatars-TotNhalos,avatardup);
  return avatardup;
}

MyIDtype assignAvatar()
{
  MyIDtype i,j,loop1,loop2,avatardup;
  avatardup = 0;
  for(i=0;i<output.nHalos;i++)
    {
      for(j=0;j<output.progs[i].nProgs;j++)
	{
	  if(HaloTable[output.progs[i].progID[j]].ProgAvatarFlag == 0)
	    {
	      HaloTable[output.progs[i].progID[j]].ProgAvatarFlag = 1;
	    }
	  else if(HaloTable[output.progs[i].progID[j]].ProgAvatarFlag == 1)
	    {
	      //printf("duplicate prog: %llu  %llu\n",IDmap[output.progs[i].progID[j]],IDmap[output.progs[i].progID[j]]);
	      HaloTable[output.progs[i].progID[j]].nAvatars++;
	      TotNavatars++;
	      HaloTable[output.progs[i].progID[j]].AvatarList = realloc(HaloTable[output.progs[i].progID[j]].AvatarList, HaloTable[output.progs[i].progID[j]].nAvatars*sizeof(MyIDtype));
	      Avatar = realloc(Avatar, TotNavatars*sizeof(MyIDtype));
	      HaloTable[output.progs[i].progID[j]].AvatarList[HaloTable[output.progs[i].progID[j]].nAvatars-1]=TotNavatars-1;
	      Avatar[TotNavatars-1] = HaloTable[output.progs[i].progID[j]].AvatarList[0];
	      //printf("duplicate prog: %llu  %llu\n",IDmap[output.progs[i].progID[j]],IDmap[Avatar[TotNavatars-1]]);
	    }
	  else
	    {
	      exit(0);
	    }
	}
    }
  for(i=0;i<TotNhalos;i++)
    {
      if(HaloTable[i].nAvatars > 1)
      	{
	  avatardup++;
      // printf("%f\t%llu\n",HaloTable[i].Mvir,HaloTable[i].nAvatars);
	}
    }
  printf("Total: TotNhalos = %llu, TotNavatar = %llu,  ExtraAvatars = %llu, MultipleAvatars = %llu\n",TotNhalos,TotNavatars,TotNavatars-TotNhalos,avatardup);
  return avatardup;
}

