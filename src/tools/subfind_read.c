#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <errno.h>
#define GENS 2
#define MAXSTRING 1024

#define MINCOMMON      10   // we only cross-correlate haloes if they at least share MINCOMMON particles
#define NSATinHOSTfrac 0.5  // at least NSATinHOSTfrac the subhaloes particles need to be inside the prospective host

//#define SO_VEL_DISPERSIONS

//#define DEBUG
//#define MTREE_SELF

// writes output that readily allows to find mergers (note, the _idx file will be replaced by _merger)
//#define MERGER_RATIO   0.25

#define CLUES_WM3

typedef struct HALOS *HALOptr;
typedef struct HALOS
{
  long   npart;
  long  *Pid;
}HALOS;

typedef struct PARTS *PARTptr;
typedef struct PARTS
{
  long  nhalos;
  long *Hid;
}PARTS;

typedef struct MTREE
{
  long unsigned id[2];
  long unsigned npart[2];
  long unsigned common;
} MTREE;

typedef unsigned int MyIDType;

struct halo_catalogue
{
  long long TotNids;
  MyIDType *SubOffset;
  MyIDType *IdList;
  int TotNsubhalos;
  int TotNgroups;
  int *SubLen;
  int *SubParentHalo;
  int *IdToHalo;
  int *GroupNsubs;
  int *GroupFirstSub;
  float *SubhaloMass;
  float *SubhaloPos;
  float *SubhaloVel;
  float *SubhaloCM;
  float *SubhaloSpin;
  float *SubhaloVelDisp;
  float *SubhaloVmax;
  float *SubhaloVmaxRad;
  float *SubhaloHalfMass;
  MyIDType *SubhaloMostBoundID;
  int *SubhaloGrNr;  
  struct descendant_data
  {
    int SnapNum[GENS];
    int HaloIndex[GENS];
#ifdef SKIP_BY_WEIGHT
    float Weight[GENS];
#endif
  }    
    *Descendant;
  int *CountProgenitors;
} CatA;
  

HALOptr halos[2];
PARTptr parts[2];
long    *main_progeniter;
long    nHalos[2];
long    PidMax=-1;
long    PidMin=1234567890;
double *merging_rate, *delta_out, *forming_rate;
MTREE *mtree;   
long unsigned  nlines;       // number of lines in *_mtree file
long unsigned  nstree;       // for how many halos to write *_stree file
long          *ihost;        // assigns one host to each halo

  
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void load_subhalo_catalogue(int num, struct halo_catalogue *cat);
char  OutputDir[1000];
void *mymalloc(size_t n);
void myfree(void *ptr);
void long_to_str(char *s, long long n);
void reassign_ids(MyIDType N, MyIDType * ids);

long total_halos(char filename[MAXSTRING]);
int lost_particles(long prog_halo, long merged_halo, int step);
int read_particles(char filename[MAXSTRING], int isimu);
int particle_halo_mapping(int isimu);

int main(int argc, char **argv)
{
  long long i,j,k,l;
  sprintf(OutputDir,"/export/virgo/Boyd/Millgas/62.5_SF_0.2/");
  load_subhalo_catalogue(61, &CatA);
  printf("%d\n",CatA.TotNgroups);
  
  for(i=0;i<CatA.TotNsubhalos;i++)
    {
      printf("%d => %f %f %f => %d\n",i,CatA.SubhaloPos[3*i]*1000.,CatA.SubhaloPos[3*i+1]*1000.,1000*CatA.SubhaloPos[3*i+2],CatA.SubhaloGrNr[i]);
    }
  /*
  halos[0] = (HALOptr) calloc(1, sizeof(HALOS));
  read_particles("/export/virgo/Boyd/Millgas/62.5_AHF_200/62.5_dm_.z0.000.AHF_particles", 0);

  parts[0] = (PARTptr) calloc(PidMax+1, sizeof(PARTS));
  particle_halo_mapping(0);
  printf("groupd %d has %d subhalos , %d as first sub\n",i,CatA.GroupNsubs[i],CatA.GroupFirstSub[i]);
  for(j=CatA.GroupFirstSub[i];j<CatA.GroupFirstSub[i]+CatA.GroupNsubs[i];j++)
    {
      printf("halo %d has %d particles: subparenthalo %d\n",j, CatA.SubLen[j],CatA.SubParentHalo[i]);
      for(k=CatA.SubOffset[j]; k< CatA.SubLen[j]+CatA.SubOffset[j];k++ )
	{
	  //  printf("\tid: %d\n",CatA.IdList[k]);
	  printf("\tid %d:%d linked to\t",CatA.IdList[k],j);
	  for(l=0;l<parts[0][CatA.IdList[k]].nhalos;l++)
	    printf("%d\t",parts[0][CatA.IdList[k]].Hid[l]);
	  printf("\n");
	}
    }
      //  }
  

  free(halos[0]);
  free(parts[0]);
  */
}


size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;
  
  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      if(feof(stream))
	printf("I/O error (fread) has occured: end of file\n");
      else
	printf("I/O error (fread) has occured: %s\n", strerror(errno));
      fflush(stdout);
      exit(778);
    }
  return nread;
} 


void load_subhalo_catalogue(int num, struct halo_catalogue *cat)
{
  int i,iboyd, ngroups, nids, nFiles, nsubhalos, subcount, groupcount,offset;
  MyIDType idcount;
  char buf[1000];
  FILE *fd;
  unsigned int j, *tmp;

  nFiles = 1;
  subcount = 0;
  idcount = 0;
  groupcount = 0;

  for(i = 0; i < nFiles; i++)
    {
      sprintf(buf, "%s/groups_%03d/subhalo_tab_%03d.%d", OutputDir, num, num, i);
      printf( "%s/groups_%03d/subhalo_tab_%03d.%d\n", OutputDir, num, num, i);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
      if(i == 0 || i == nFiles - 1)
	printf("reading '%s'\n", buf);
      if(i == 1)
	printf("...to...\n");


      my_fread(&ngroups, sizeof(int), 1, fd);
      printf("ngroup = %d\n",ngroups);
      my_fread(&cat->TotNgroups, sizeof(int), 1, fd);
      printf("TotNgroup = %d\n",cat->TotNgroups);
      my_fread(&nids, sizeof(int), 1, fd);
      printf("nids = %d\n",nids);
      my_fread(&cat->TotNids, sizeof(long long), 1, fd);
      printf("TotNids = %d\n",cat->TotNids);
      my_fread(&nFiles, sizeof(int), 1, fd);
      my_fread(&nsubhalos, sizeof(int), 1, fd);
      printf("nsubhalos = %d\n",nsubhalos);
      my_fread(&cat->TotNsubhalos, sizeof(int), 1, fd);
      printf("TotNsubhalos = %d\n",cat->TotNsubhalos);
  

      if(i == 0)
	{
	  cat->IdList = mymalloc(sizeof(MyIDType) * cat->TotNids);
	  cat->SubLen = mymalloc(sizeof(int) * cat->TotNsubhalos);
	  cat->SubOffset = mymalloc(sizeof(MyIDType) * cat->TotNsubhalos);

	  cat->SubParentHalo = mymalloc(sizeof(int) * cat->TotNsubhalos);
	  cat->SubhaloMass = mymalloc(sizeof(float) * cat->TotNsubhalos);

	  cat->SubhaloPos = mymalloc(3 * sizeof(float ) * cat->TotNsubhalos);
	  cat->SubhaloVel = mymalloc(3 * sizeof(float ) * cat->TotNsubhalos);
	  cat->SubhaloCM = mymalloc(3 * sizeof(float ) * cat->TotNsubhalos);
	  cat->SubhaloSpin = mymalloc(3 * sizeof(float ) * cat->TotNsubhalos);


	  cat->SubhaloVelDisp = mymalloc(sizeof(float) * cat->TotNsubhalos);
	  cat->SubhaloVmax = mymalloc(sizeof(float) * cat->TotNsubhalos);
	  cat->SubhaloVmaxRad = mymalloc(sizeof(float) * cat->TotNsubhalos);
	  cat->SubhaloHalfMass = mymalloc(sizeof(float) * cat->TotNsubhalos);

	  cat->SubhaloMostBoundID = mymalloc(sizeof(MyIDType) * cat->TotNsubhalos);
	  cat->SubhaloGrNr = mymalloc(sizeof(int) * cat->TotNsubhalos);

	  cat->GroupNsubs = mymalloc(sizeof(int) * cat->TotNgroups);
	  cat->GroupFirstSub = mymalloc(sizeof(int) * cat->TotNgroups);
	  cat->Descendant = mymalloc(sizeof(struct descendant_data) * cat->TotNsubhalos);
	  cat->CountProgenitors = mymalloc(sizeof(int) * cat->TotNsubhalos);
	}

      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupLen  */
      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupOffset  */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  GroupMass  */
      fseek(fd, 3 * sizeof(float) * ngroups, SEEK_CUR);	/* skip  GroupPos */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_M_Mean200 */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_R_Mean200 */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_M_Crit200 */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_R_Crit200 */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_M_TopHat200 */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_R_TopHat200 */
#ifdef SO_VEL_DISPERSIONS
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_VelDisp_Mean200 */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_VelDisp_Crit200 */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_VelDisp_TopHat200 */
#endif
      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupContaminationCount */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  GroupContaminationMass */
      //fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupNsubs */
      my_fread(&cat->GroupNsubs[groupcount], sizeof(int), ngroups, fd);
      //fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupFirstsub */
      my_fread(&cat->GroupFirstSub[groupcount], sizeof(int), ngroups, fd);
      
      my_fread(&cat->SubLen[subcount], sizeof(int), nsubhalos, fd);

      tmp = mymalloc(sizeof(int) * nsubhalos);
      my_fread(tmp, sizeof(int), nsubhalos, fd);
      for(j = 0; j < nsubhalos; j++)
	cat->SubOffset[subcount + j] = tmp[j];	/* copy it to 64 bit if needed */

      myfree(tmp);
      //myfree(buf);


      my_fread(&cat->SubParentHalo[subcount], sizeof(int), nsubhalos, fd);
      my_fread(&cat->SubhaloMass[subcount], sizeof(float), nsubhalos, fd);
      //fseek(fd, sizeof(float) * nsubhalos, SEEK_CUR);	/* skip  SubhaloMass */

      my_fread(&cat->SubhaloPos[3 * subcount], 3 * sizeof(float), nsubhalos, fd);
      my_fread(&cat->SubhaloVel[3 * subcount], 3 * sizeof(float), nsubhalos, fd);

      //fseek(fd, 3 * sizeof(float) * nsubhalos, SEEK_CUR);	/* skip  SubhaloCM */
      my_fread(&cat->SubhaloCM[3 * subcount], 3 * sizeof(float), nsubhalos, fd);

      my_fread(&cat->SubhaloVelDisp[subcount], sizeof(float), nsubhalos, fd);
      my_fread(&cat->SubhaloVmax[subcount], sizeof(float), nsubhalos, fd);

      my_fread(&cat->SubhaloVmaxRad[subcount], sizeof(float), nsubhalos, fd);

      //fseek(fd, sizeof(float) * nsubhalos, SEEK_CUR);	/* skip  SubhaloVmaxRad */

      my_fread(&cat->SubhaloHalfMass[subcount], sizeof(float), nsubhalos, fd);

      
      my_fread(&cat->SubhaloMostBoundID[subcount], sizeof(MyIDType), nsubhalos, fd);

      //fseek(fd, sizeof(int) * nsubhalos, SEEK_CUR);	/* skip  GrNr */
      my_fread(&cat->SubhaloGrNr[subcount], sizeof(int), nsubhalos, fd);
      
      fclose(fd);

      subcount += nsubhalos;
      groupcount += ngroups;
    }

  /* now check whether we had a 32-bit overflow in subhalo offset, and fix it if needed */

  if(sizeof(MyIDType) > 4)
    {
      MyIDType previous_offset = 0, add = 0;

      for(i = 0; i < cat->TotNsubhalos; i++)
	{
	  if(cat->SubOffset[i] < previous_offset)
	    {
	      printf("deteced 32-bit overflow, correcting it\n");
	      add += (((long long) 1) << 32);
	    }
	  previous_offset = cat->SubOffset[i];

	  cat->SubOffset[i] += add;
	}
    }


  long_to_str(buf, cat->TotNids);

  printf("cat->TotNsubhalos = %d\n", cat->TotNsubhalos);
  printf("cat->TotNids = %s\n", buf);


  for(i = 0; i < nFiles; i++)
    {
      sprintf(buf, "%s/groups_%03d/subhalo_ids_%03d.%d", OutputDir, num, num, i);
      printf("%s/groups_%03d/subhalo_ids_%03d.%d\n", OutputDir, num, num, i);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
      if(i == 0 || i == nFiles - 1)
	printf("reading '%s'\n", buf);
      if(i == 1)
	printf("...to...\n");

      my_fread(&ngroups, sizeof(int), 1, fd);
      my_fread(&cat->TotNgroups, sizeof(int), 1, fd);
      my_fread(&nids, sizeof(int), 1, fd);
      my_fread(&cat->TotNids, sizeof(long long), 1, fd);
      my_fread(&nFiles, sizeof(int), 1, fd);
      my_fread(&offset, sizeof(int), 1, fd);

      printf("read all\n");
      my_fread(&cat->IdList[idcount], sizeof(MyIDType), nids, fd);
      fclose(fd);

      idcount += nids;
    }


  //reassign_ids(cat->TotNids, cat->IdList);
}

void *mymalloc(size_t n)
{
  void *p;

  if(!(p = malloc(n)))
    {
      if(n)
	{
	  printf("Failed to allocate memory for %u bytes.\n", (int) n);
	  exit(2);
	}
    }

  return p;
}


void myfree(void *ptr)
{
  free(ptr);
}

void reassign_ids(MyIDType N, MyIDType * ids)
{
#ifdef IDS_HAVE_GAPS

  long long i, j, offset, NN;
  int tid, nthreads;
  struct twoids *TwoIDs;


  printf("reassign IDs...\n");
  fflush(stdout);

#if defined(_OPENMP)
#pragma omp parallel private(tid, nthreads, offset, NN, i, j, TwoIDs) shared(IdSnapTable)
#endif
  {
#if defined(_OPENMP)
    tid = omp_get_thread_num();
    nthreads = omp_get_max_threads();
    
    offset = tid * (N / nthreads);
    NN = (N / nthreads);
    
    if(nthreads > 1 && tid == (nthreads - 1))
      {
	NN = N - offset;
      }
    printf("tid=%d offset=%lld NN=%lld\n", tid, offset, NN); 

#else
    NN = N;
    offset = 0;
    tid = 0;
    nthreads = 1;
#endif

    TwoIDs = mymalloc(NN * sizeof(struct twoids));


    for(i = 0; i < NN; i++)
      {
	TwoIDs[i].id = ids[i + offset];
	TwoIDs[i].ord = i;
      }

    qsort(TwoIDs, NN, sizeof(struct twoids), sort_twoids_id);

    /* now assign */

    j = 0;
    for(i = 0; i < NN; i++)
      {
	while(IdSnapTable[j] < TwoIDs[i].id && j < (TotNumPart - 1))
	  j++;

	if(IdSnapTable[j] != TwoIDs[i].id)
	  {
	    printf("ID mismatch found?\n");
	    exit(1);
	  }

	TwoIDs[i].id = j;
      }


    /* sort back */
    qsort(TwoIDs, NN, sizeof(struct twoids), sort_twoids_ord);


    for(i = 0; i < NN; i++)
      {
	ids[i + offset] = TwoIDs[i].id;
      }
    myfree(TwoIDs);
  }

  printf("done\n");
  fflush(stdout);

#else
  
  signed long long i;
 
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for(i=0; i< N; i++)
    ids[i] -= 1;

#endif

}

void long_to_str(char *s, long long n)
{
  if(n >= 1000000000)
    sprintf(s, "%d%09d", (int) (n / 1000000000), (int) (n % 1000000000));
  else
    sprintf(s, "%d", (int) n);
}
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



#ifdef Boyd

void load_subhalo_catalogue(int num)
{
  int i, ngroups, nids, nFiles, nsubhalos, subcount;
  int groupcount, filenr, ncount;
  int subgr, gr, nh, sc, gr_nh;
  char buf[1000];
  FILE *fd;
  int *nsubPerHalo, *subLen, *descendant_haloindex, *descendant_snapnum, *filenrOfHalo, *subhaloindex;
  float *halo_M_Mean200, *halo_M_Crit200, *halo_M_TopHat;
  float *subpos, *subvel, *subveldisp, *subvmax, *subspin, *subhalfmass;
  MyIDType *subMostBoundID;

#ifdef SAVE_MASS_TAB
  float *submasstab;
#endif

  printf("Catalogue num=%d\n", num);


  nsubPerHalo = mymalloc(sizeof(int) * Cats[num].TotNgroups);
  subLen = mymalloc(sizeof(int) * Cats[num].TotNsubhalos);
  descendant_haloindex = mymalloc(sizeof(int) * Cats[num].TotNsubhalos);
  descendant_snapnum = mymalloc(sizeof(int) * Cats[num].TotNsubhalos);
  filenrOfHalo = mymalloc(sizeof(int) * Cats[num].TotNsubhalos);
  subhaloindex = mymalloc(sizeof(int) * Cats[num].TotNsubhalos);
  
  halo_M_Mean200 = mymalloc(sizeof(float) * Cats[num].TotNgroups);
  halo_M_Crit200 = mymalloc(sizeof(float) * Cats[num].TotNgroups);
  halo_M_TopHat = mymalloc(sizeof(float) * Cats[num].TotNgroups);
  
  subpos = mymalloc(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subvel = mymalloc(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subveldisp = mymalloc(sizeof(float) * Cats[num].TotNsubhalos);
  subvmax = mymalloc(sizeof(float) * Cats[num].TotNsubhalos);
  subspin = mymalloc(3 * sizeof(float) * Cats[num].TotNsubhalos);
  subMostBoundID = mymalloc(sizeof(MyIDType) * Cats[num].TotNsubhalos);
  subhalfmass = mymalloc(sizeof(float) * Cats[num].TotNsubhalos);
#ifdef SAVE_MASS_TAB
  submasstab = mymalloc(6 * sizeof(float) * Cats[num].TotNsubhalos);
#endif


  subcount = 0;
  groupcount = 0;

  nFiles = 1;

  for(filenr = 0; filenr < nFiles; filenr++)
    {
      sprintf(buf, "%s/groups_%03d/subhalo_tab_%03d.%d", OutputDir, num, num, filenr);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}
      long long totNids;

      my_fread(&ngroups, sizeof(int), 1, fd);
      my_fread(&Cats[num].TotNgroups, sizeof(int), 1, fd);
      my_fread(&nids, sizeof(int), 1, fd);
      my_fread(&totNids, sizeof(long long), 1, fd);
      my_fread(&nFiles, sizeof(int), 1, fd);
      my_fread(&nsubhalos, sizeof(int), 1, fd);
      my_fread(&Cats[num].TotNsubhalos, sizeof(int), 1, fd);

#ifdef Boyd

      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupLen  */
      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupOffset  */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  GroupMass  */
      fseek(fd, 3 * sizeof(float) * ngroups, SEEK_CUR);	/* skip  GroupPos */

      my_fread(&halo_M_Mean200[groupcount], sizeof(float), ngroups, fd);
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_R_Mean200 */

      my_fread(&halo_M_Crit200[groupcount], sizeof(float), ngroups, fd);
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_R_Crit200 */

      my_fread(&halo_M_TopHat[groupcount], sizeof(float), ngroups, fd);
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_R_TopHat200 */

#ifdef FLAG_GROUP_VELDISP
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_VelDisp_Mean200 */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_VelDisp_Crit200 */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_VelDisp_TopHat200 */
#endif
      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupContaminationCount */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  GroupContaminationMass */




      my_fread(&nsubPerHalo[groupcount], sizeof(int), ngroups, fd);

      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupFirstsub */

      my_fread(&subLen[subcount], sizeof(int), nsubhalos, fd);


      fseek(fd, sizeof(int) * nsubhalos, SEEK_CUR);	/* skip  SubOffset */
      fseek(fd, sizeof(int) * nsubhalos, SEEK_CUR);	/* skip  SubParenthalo */

      fseek(fd, sizeof(float) * nsubhalos, SEEK_CUR);	/* skip  SubhaloMass */

      my_fread(&subpos[3 * subcount], 3 * sizeof(float), nsubhalos, fd);
      my_fread(&subvel[3 * subcount], 3 * sizeof(float), nsubhalos, fd);

      fseek(fd, 3 * sizeof(float) * nsubhalos, SEEK_CUR);	/* skip  SubhaloCM */

      my_fread(&subspin[3 * subcount], 3 * sizeof(float), nsubhalos, fd);
      my_fread(&subveldisp[subcount], sizeof(float), nsubhalos, fd);
      my_fread(&subvmax[subcount], sizeof(float), nsubhalos, fd);

      fseek(fd, sizeof(float) * nsubhalos, SEEK_CUR);	/* skip  SubhaloVmaxRad */
      my_fread(&subhalfmass[subcount], sizeof(float), nsubhalos, fd);


      my_fread(&subMostBoundID[subcount], sizeof(MyIDType), nsubhalos, fd);

      fseek(fd, sizeof(int) * nsubhalos, SEEK_CUR);	/* skip  GrNr */

#ifdef SAVE_MASS_TAB
      my_fread(&submasstab[6 * subcount], 6 * sizeof(float), nsubhalos, fd);
#endif
      fclose(fd);


      for(subgr = 0; subgr < nsubhalos; subgr++)
	filenrOfHalo[subcount + subgr] = filenr;

      for(subgr = 0; subgr < nsubhalos; subgr++)
	subhaloindex[subcount + subgr] = subgr;

      subcount += nsubhalos;
      groupcount += ngroups;

#endif
    }


  if(num < LastSnapShotNr)
    {
      sprintf(buf, "%s/treedata/sub_desc_%03d", OutputDir, num);
      //sprintf(buf, "%s/treedata/sub_desc_sf%d_%03d", OutputDir, SnapSkipFac, num);
      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s'\n", buf);
	  exit(1);
	}

      my_fread(&ncount, sizeof(int), 1, fd);
      my_fread(descendant_haloindex, sizeof(int), Cats[num].TotNsubhalos, fd);
      my_fread(descendant_snapnum, sizeof(int), Cats[num].TotNsubhalos, fd);

      fclose(fd);
    }

  nh = FirstHaloInSnap[num];
  sc = 0;

  for(gr = 0; gr < Cats[num].TotNgroups; gr++)
    {
      for(subgr = 0, gr_nh = nh; subgr < nsubPerHalo[gr]; subgr++, sc++, nh++)
	{
	  Halo[nh].FirstHaloInFOFgroup = gr_nh;
	  if(subgr == nsubPerHalo[gr] - 1)
	    Halo[nh].NextHaloInFOFgroup = -1;
	  else
	    Halo[nh].NextHaloInFOFgroup = nh + 1;

	  if(num < LastSnapShotNr)
	    {
	      if(descendant_haloindex[sc] >= 0)
		Halo[nh].Descendant = FirstHaloInSnap[descendant_snapnum[sc]] + descendant_haloindex[sc];
	      else
		Halo[nh].Descendant = -1;
	    }
	  else
	    Halo[nh].Descendant = -1;

	  Halo[nh].FirstProgenitor = -1;
	  Halo[nh].NextProgenitor = -1;

	  /* assign properties */

	  Halo[nh].Len = subLen[sc];

	  if(subgr == 0)
	    {
	      Halo[nh].M_Mean200 = halo_M_Mean200[gr];
	      Halo[nh].M_Crit200 = halo_M_Crit200[gr];
	      Halo[nh].M_TopHat = halo_M_TopHat[gr];
	    }
	  else
	    {
	      Halo[nh].M_Mean200 = 0;
	      Halo[nh].M_Crit200 = 0;
	      Halo[nh].M_TopHat = 0;
	    }


	  for(i = 0; i < 3; i++)
	    {
	      Halo[nh].Pos[i] = subpos[3 * sc + i];
	      Halo[nh].Vel[i] = subvel[3 * sc + i];
	      Halo[nh].Spin[i] = subspin[3 * sc + i];
	    }
	  Halo[nh].VelDisp = subveldisp[sc];
	  Halo[nh].Vmax = subvmax[sc];
	  Halo[nh].MostBoundID = subMostBoundID[sc];


	  /* store position of halo in subfind output */

	  Halo[nh].SnapNum = num;
	  Halo[nh].FileNr = filenrOfHalo[sc];
	  Halo[nh].SubhaloIndex = subhaloindex[sc];
	  Halo[nh].SubhalfMass = subhalfmass[sc];

	  /* auxiliary stuff */

	  HaloAux[nh].UsedFlag = 0;

#ifdef SAVE_MASS_TAB
	  for(i = 0; i < 6; i++)
	    Halo[nh].SubMassTab[i] = submasstab[6 * sc + i];
#endif
	}
    }


  for(gr = 0; gr < nh; gr++)
    {
      if(Halo[gr].NextHaloInFOFgroup == gr)
	{
	  printf("bummer! %d\n", gr);
	}
    }




#ifdef SAVE_MASS_TAB
  myfree(submasstab);
#endif
  myfree(subhalfmass);
  myfree(subMostBoundID);
  myfree(subspin);
  myfree(subvmax);
  myfree(subveldisp);
  myfree(subvel);
  myfree(subpos);

  myfree(halo_M_TopHat);
  myfree(halo_M_Crit200);
  myfree(halo_M_Mean200);

  myfree(subhaloindex);
  myfree(filenrOfHalo);
  myfree(descendant_snapnum);
  myfree(descendant_haloindex);
  myfree(subLen);
  myfree(nsubPerHalo);
}

#endif
