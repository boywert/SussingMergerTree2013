#include "allvars.h"

void SubtreeFilter();
void resetHaloTable();

char FilePrefix[MAXSTRING];
char HBTFilePrefix[MAXSTRING];
char FolderName[MAXSTRING];
char HBTFolderName[MAXSTRING];
char TreeFolder[MAXSTRING];
char ppFolder[MAXSTRING];
char inputFile[MAXSTRING];
char *trimwhitespace(char *str);
char gadgetfolder[MAXSTRING];
char gadgetPrefix[MAXSTRING];
char filteredFolder[MAXSTRING];
int main(int argc, char **argv)
{
  int isTree;
  MyIDtype i;
  char buffer[MAXSTRING];
  FILE* fp;
  isTree = 1;
  //FilePrefix = calloc(MAXSTRING,sizeof(char));
  //FolderName = calloc(MAXSTRING,sizeof(char));
  //  strcpy(FolderName,"/export/virgo/Boyd/Millgas/62.5_ALEX");
  //  strcpy(FilePrefix,"62.5_dm_");
  fp = fopen("config.txt","r");

  fgets(buffer,MAXSTRING,fp);
  sprintf(FolderName,"%s",trimwhitespace(buffer));

  fgets(buffer,MAXSTRING,fp);
  sprintf(FilePrefix,"%s",trimwhitespace(buffer));

  fgets(buffer,MAXSTRING,fp);
  sprintf(SnapTimeFile,"%s",trimwhitespace(buffer));

  fgets(buffer,MAXSTRING,fp);
  sprintf(TreeFolder,"%s",trimwhitespace(buffer));
  
  fgets(buffer,MAXSTRING,fp);
  
  fgets(buffer,MAXSTRING,fp);
  sprintf(ppFolder,"%s",trimwhitespace(buffer));
  fclose(fp);
    
  //  sprintf(SnapTimeFile,"%s/%s",FolderName,"data_snaplist.txt");
  getSnapTime();

  if(argc < 2)
    {
      printf("\n  usage: %s <outputfile>\n",argv[0]);
      printf("  <outputfile>        output filename\n\n");
      exit(1);
    }

  if(argc == 4)
    {
      strcpy(delFile,argv[2]);
      strcpy(addFile,argv[3]);
    }
  if(argc == 3)
    {
      strcpy(addFile,argv[2]);
    }
  strcpy(inputFile,argv[1]);
  sprintf(HBTFolderName,"%s",addFile);
  sprintf(HBTFilePrefix,"%s",FilePrefix);
  if(readRawOutputs(argv[1]) == 0)
    {
      fprintf(stderr,"The output: %s is valid\n",argv[1]);
    }
  else
    {
      fprintf(stderr,"The output: %s is invalid\n",argv[1]);
    }
  qsort_r(output.progs,output.nHalos, sizeof(struct progtable), compareHaloID);

  //SubtreeFilter();
  if(assignAvatarMain() != 0)
    {
      printf("Not trees\n");
      isTree = 0;
    }
  qsort(output.progs,output.nHalos, sizeof(struct progtable), compareHaloID);
#ifdef datasetII
  for(i=0;i<TotNhalos;i++)
    {
      if(SubTree[i] < NULLPOINT)
	HaloTable[i].Mvir = 0.;
    }
#endif
  (void) printouttrees();
  //printoutparticles_binary();
  //exit(0);
  //mean_main_branch_mass();
  main_branch_analysis(0.,1e11,100);  
  main_branch_analysis(2e11,5e11,500);
  main_branch_analysis(1e12,1e30,1000);
  if(isTree)
    {
      merger_analysis(0.,1e11,100);
      //merger_analysis(2e11,5e11,500);
      //merger_analysis(1e12,1e30,1000);
    }
  //snapshot_stats(61,0.,1e11,100);
  //snapshot_stats(61,2e11,5e11,500);
  //snapshot_stats(61,1e12,1e30,1000);
  //calculateMassFunction();
  return 0;
}

char *trimwhitespace(char *str)
{
  char *end;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
    return str;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  // Write new null terminator
  *(end+1) = 0;

  return str;
}

/*------------------------------------------------------------------------------
Utilities 
-------------------------------------------------------------------------------*/
/*
void SubtreeFilter()
{
  MyIDtype i,j,k,l,MainProg,current,structuredup;
  printf("Use Subtree Filter ...\n");
  structuredup = 0;
  for(i=0;i<output.nHalos;i++)
    {
      //printf("i = %llu\n",i);
      if(output.progs[i].nProgs > 1)
	{
	  for(j=0;j<output.progs[i].nProgs;j++)
	    {
	      //printf("i = %llu, j = %llu\n",i,j);
	      current = SubTree[output.progs[i].progID[j]];
	      while(current != NULLPOINT)
		{
		  //printf("i = %llu, j = %llu, current = %llu\n",i,j,current);
		  for(k=j+1;k<output.progs[i].nProgs;k++)
		    {
		      if(current == output.progs[i].progID[k])
			{
			  //fprintf(stderr,"got substructure problem\n");
			  for(l=k+1;l<output.progs[i].nProgs-1;l++)
			    {
			      output.progs[i].progID[l] = output.progs[i].progID[l+1];
			    }
			  output.progs[i].nProgs--;
			  output.progs[i].progID = realloc( output.progs[i].progID,output.progs[i].nProgs*sizeof(MyIDtype));
			  structuredup++;
			}
		    }
		  current = SubTree[current];
		}
	      for(k=j+1;k<output.progs[i].nProgs;k++)
		{
		  current = SubTree[output.progs[i].progID[k]];
		  while(current != NULLPOINT)
		    {
		      if(current == output.progs[i].progID[j])
			{
			  //fprintf(stderr,"got substructure problem\n");
			  for(l=k+1;l<output.progs[i].nProgs-1;l++)
			    {
			      output.progs[i].progID[l] = output.progs[i].progID[l+1];
			    }
			  output.progs[i].nProgs--;
			  output.progs[i].progID = realloc( output.progs[i].progID,output.progs[i].nProgs*sizeof(MyIDtype));
			  structuredup++;
			}
		      current = SubTree[current];
		    }
		}
	    }
	}
    }
}





void resetHaloTable()
{
  MyIDtype i;
  for(i=0;i<TotNhalos;i++)
    {
      HaloTable[i].ProgAvatarFlag = 0;
    }
}


*/
