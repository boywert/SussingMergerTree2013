#include "allvars.h"


char delFile[MAXSTRING];
char addFile[MAXSTRING];
char SnapTimeFile[MAXSTRING];
MyIDtype SnapNhalos[NSNAPS];
float snapTime[NSNAPS];
float snapTimeYear[NSNAPS];
float expansion_factor[NSNAPS];
MyIDtype TotNhalos;
MyIDtype TotNhalosUsed;
MyIDtype *IDmap;
MyIDtype *SubTree;


struct SNAPSHOT_STATS snap_stats[NSNAPS];
struct HALOPROPS *HaloTable;
struct HALOPROPS *HBThaloTable;
struct HBT_halos *HBT;


void read_particles(unsigned int slotid);
void read_particles_binary();

void getSnapTime()
{
  int i,dummyint;
  FILE* fp;
  char line[MAXSTRING];
  float dummyfloat;
  printf("Start getting snapshot time: %s\n",SnapTimeFile);
  fp = fopen(SnapTimeFile, "r");
  fgets(line,MAXSTRING,fp);
  for(i=0;i<NSNAPS;i++)
    {
      fscanf(fp,"%d %g %g %g %g",
	     &dummyint,
	     &(expansion_factor[i]),
	     &dummyfloat,
	     &(snapTime[i]),
	     &(snapTimeYear[i])
	     );
    }
  
}

void makeIDmap()
{
  int len;
  unsigned int iFile, i;
  MyIDtype currentHalo, nHalos,iHalo,hosthalo;
  char filename[MAXSTRING];
  char buffer[MAXSTRING];
  char dummystr[MAXSTRING];
  FILE* fp;
  TotNhalos = 0;
  currentHalo = 0;
  for(iFile=FIRSTSNAP;iFile<=LASTSNAP;iFile++)
    {
      //printf("iFile = %d\n",iFile);
      strncpy(filename,"",sizeof(filename));
      (void) getFilename(filename,iFile);
      //printf("%s\n",strcat(filename,".AHF_particles"));
      fp = fopen(strcat(filename,".AHF_particles"), "r");
      fscanf(fp,"%llu",&(SnapNhalos[iFile]));
      fclose(fp);
      TotNhalos += SnapNhalos[iFile];
      //printf("%d snap = %llu halos\n",iFile,SnapNhalos[iFile]);
    }
  if(TotNhalos <= MAXUSEABLE)
    {
      IDmap = calloc(TotNhalos,sizeof(MyIDtype));
      HaloTable = calloc(TotNhalos,sizeof(struct HALOPROPS));
      TotNavatars = TotNhalos;
      TotNhalosUsed = TotNhalos;
      Avatar= calloc(TotNhalos,sizeof(MyIDtype));
    }
  else
    {
      fprintf(stderr,"Error: Total halo exceed the maximum:%llu\nExiting...\n",MAXUSEABLE);
      exit(0);
    }
  for(iFile=FIRSTSNAP;iFile<=LASTSNAP;iFile++)
    {
      strncpy(filename,"",sizeof(filename));
      (void) getFilename(filename,iFile);
      fp = fopen(strcat(filename,".AHF_halos"), "r");
      fgets(buffer,MAXSTRING,fp);
      for(iHalo=0;iHalo<SnapNhalos[iFile];iHalo++)
	{
	  fscanf(fp,"%llu %llu %llu %g %llu %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g",
		 &(HaloTable[currentHalo].ID),
		 &(HaloTable[currentHalo].hostHalo),
		 &(HaloTable[currentHalo].numSubStruct),
		 &(HaloTable[currentHalo].Mvir),
		 &(HaloTable[currentHalo].npart),
		 &(HaloTable[currentHalo].Xc),
		 &(HaloTable[currentHalo].Yc),
		 &(HaloTable[currentHalo].Zc),
		 &(HaloTable[currentHalo].VXc),
		 &(HaloTable[currentHalo].VYc),
		 &(HaloTable[currentHalo].VZc),
		 &(HaloTable[currentHalo].Rvir), 
		 &(HaloTable[currentHalo].Rmax),
		 &(HaloTable[currentHalo].r2),
		 &(HaloTable[currentHalo].mbp_offset),
		 &(HaloTable[currentHalo].com_offset),
		 &(HaloTable[currentHalo].Vmax),
		 &(HaloTable[currentHalo].v_esc),
		 &(HaloTable[currentHalo].sigV),
		 &(HaloTable[currentHalo].lambda),
		 &(HaloTable[currentHalo].lambdaE),
		 &(HaloTable[currentHalo].Lx),
		 &(HaloTable[currentHalo].Ly),
		 &(HaloTable[currentHalo].Lx),
		 &(HaloTable[currentHalo].b),
		 &(HaloTable[currentHalo].c),
		 &(HaloTable[currentHalo].Eax),
		 &(HaloTable[currentHalo].Eay),
		 &(HaloTable[currentHalo].Eaz),
		 &(HaloTable[currentHalo].Ebx), 
		 &(HaloTable[currentHalo].Eby),
		 &(HaloTable[currentHalo].Ebz),
		 &(HaloTable[currentHalo].Ecx), 
		 &(HaloTable[currentHalo].Ecy), 
		 &(HaloTable[currentHalo].Ecy),
		 &(HaloTable[currentHalo].ovdens),
		 &(HaloTable[currentHalo].nbins),
		 &(HaloTable[currentHalo].fMhires),
		 &(HaloTable[currentHalo].Ekin),
		 &(HaloTable[currentHalo].Epot),
		 &(HaloTable[currentHalo].SurfP),
		 &(HaloTable[currentHalo].Phi0),
		 &(HaloTable[currentHalo].cNFW) 
		 );
	  IDmap[currentHalo] = HaloTable[currentHalo].ID;
	  HaloTable[currentHalo].AHFID = HaloTable[currentHalo].ID; 
	  HaloTable[currentHalo].SnapID = iFile;
	  HaloTable[currentHalo].nAvatars = 1;
	  HaloTable[currentHalo].AvatarList = calloc(HaloTable[currentHalo].nAvatars,sizeof(MyIDtype));
	  HaloTable[currentHalo].AvatarList[0] = currentHalo;
	  Avatar[currentHalo] = currentHalo;
	  HaloTable[currentHalo].ProgAvatarFlag = 0;
	  HaloTable[currentHalo].TroubleFlag = 0;

	  //printf("%llu : %llu\n", currentHalo,IDmap[currentHalo]);
	  currentHalo++;
	}
      fclose(fp);
      
    }

  //for(iFile=FIRSTSNAP;iFile<=LASTSNAP;iFile++)
  //  {
  //    read_particles(iFile);
  //  }
  if(TotNhalos != (currentHalo))
    {
      fprintf(stderr,"Error: Problem in halo catalogues: %llu,%llu\nExiting...\n",TotNhalos,currentHalo);
      exit(0);
    }
  SubTree = calloc(TotNhalos,sizeof(MyIDtype));
  for(iHalo=0;iHalo<TotNhalos;iHalo++)
    {
      SubTree[iHalo] = IDsearch(HaloTable[iHalo].hostHalo);
      //printf("halo %llu:%llu hosted in %llu: %llu\n", iHalo,IDmap[iHalo], SubTree[iHalo],HaloTable[iHalo].hostHalo);
    }
  
  for(iHalo=0;iHalo<TotNhalos;iHalo++)
    {
      hosthalo = SubTree[iHalo];
      if(hosthalo < NULLPOINT)
	{
	  if(HaloTable[hosthalo].Mvir < 1.e-30)
	    HaloTable[iHalo].Mvir = 0.;
	}
    }
  //load_particles(IDmap[1000000]);
  //load_particles(IDmap[500000]);
  if(output.outputFormat > 1.109 && output.outputFormat < 1.111)
    {
      (void) deleteHalos_v1(delFile);
      (void) addHalo_v1(addFile);
    }
  if(output.outputFormat > 1.129 && output.outputFormat < 1.131)
    {
      (void) hbtmaphalos(addFile);
#ifdef READPARTICLE
      read_particles_binary();
#endif
      /* for(iFile=FIRSTSNAP;iFile<=LASTSNAP;iFile++) */
      /* 	{ */
      /* 	  read_particles(iFile); */
      /* 	} */
      /* (void) printoutparticles_binary(); */
    }
#ifdef READPARTICLE
  if(output.outputFormat > 0.999 && output.outputFormat < 1.001)
    read_particles_binary();
#endif
  //(void) resetIDmap();
  printf("Output binary cache files\n");
  (void) get_snap_stats();
  (void) printoutprofile();

#ifdef SUBFINDOUT
  (void) makesubfindout();
#endif
}
void read_singlesnap(unsigned int snapnum)
{
  int len;
  unsigned int iFile, i;
  MyIDtype currentHalo, nHalos,iHalo,hosthalo;
  char filename[MAXSTRING];
  char buffer[MAXSTRING];
  char dummystr[MAXSTRING];
  FILE* fp;
  struct Gadget_particle *P;
  int *PIDmap;
  TotNhalos = 0;
  currentHalo = 0;
  for(iFile=snapnum;iFile<=snapnum;iFile++)
    {
      //printf("iFile = %d\n",iFile);
      strncpy(filename,"",sizeof(filename));
      (void) getFilename(filename,iFile);
      //printf("%s\n",strcat(filename,".AHF_particles"));
      fp = fopen(strcat(filename,".AHF_particles"), "r");
      fscanf(fp,"%llu",&(SnapNhalos[iFile]));
      fclose(fp);
      TotNhalos += SnapNhalos[iFile];
      //printf("%d snap = %llu halos\n",iFile,SnapNhalos[iFile]);
    }
  if(TotNhalos <= MAXUSEABLE)
    {
      IDmap = calloc(TotNhalos,sizeof(MyIDtype));
      HaloTable = calloc(TotNhalos,sizeof(struct HALOPROPS));
      TotNavatars = TotNhalos;
      TotNhalosUsed = TotNhalos;
      Avatar= calloc(TotNhalos,sizeof(MyIDtype));
    }
  else
    {
      fprintf(stderr,"Error: Total halo exceed the maximum:%llu\nExiting...\n",MAXUSEABLE);
      exit(0);
    }
  for(iFile=snapnum;iFile<=snapnum;iFile++)
    {
      strncpy(filename,"",sizeof(filename));
      (void) getFilename(filename,iFile);
      fp = fopen(strcat(filename,".AHF_halos"), "r");
      fgets(buffer,MAXSTRING,fp);
      for(iHalo=0;iHalo<SnapNhalos[iFile];iHalo++)
	{
	  //printf("iHalo = %llu\n",iHalo);
	  fscanf(fp,"%llu %llu %llu %g %llu %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g",
		 &(HaloTable[currentHalo].ID),
		 &(HaloTable[currentHalo].hostHalo),
		 &(HaloTable[currentHalo].numSubStruct),
		 &(HaloTable[currentHalo].Mvir),
		 &(HaloTable[currentHalo].npart),
		 &(HaloTable[currentHalo].Xc),
		 &(HaloTable[currentHalo].Yc),
		 &(HaloTable[currentHalo].Zc),
		 &(HaloTable[currentHalo].VXc),
		 &(HaloTable[currentHalo].VYc),
		 &(HaloTable[currentHalo].VZc),
		 &(HaloTable[currentHalo].Rvir), 
		 &(HaloTable[currentHalo].Rmax),
		 &(HaloTable[currentHalo].r2),
		 &(HaloTable[currentHalo].mbp_offset),
		 &(HaloTable[currentHalo].com_offset),
		 &(HaloTable[currentHalo].Vmax),
		 &(HaloTable[currentHalo].v_esc),
		 &(HaloTable[currentHalo].sigV),
		 &(HaloTable[currentHalo].lambda),
		 &(HaloTable[currentHalo].lambdaE),
		 &(HaloTable[currentHalo].Lx),
		 &(HaloTable[currentHalo].Ly),
		 &(HaloTable[currentHalo].Lx),
		 &(HaloTable[currentHalo].b),
		 &(HaloTable[currentHalo].c),
		 &(HaloTable[currentHalo].Eax),
		 &(HaloTable[currentHalo].Eay),
		 &(HaloTable[currentHalo].Eaz),
		 &(HaloTable[currentHalo].Ebx), 
		 &(HaloTable[currentHalo].Eby),
		 &(HaloTable[currentHalo].Ebz),
		 &(HaloTable[currentHalo].Ecx), 
		 &(HaloTable[currentHalo].Ecy), 
		 &(HaloTable[currentHalo].Ecy),
		 &(HaloTable[currentHalo].ovdens),
		 &(HaloTable[currentHalo].nbins),
		 &(HaloTable[currentHalo].fMhires),
		 &(HaloTable[currentHalo].Ekin),
		 &(HaloTable[currentHalo].Epot),
		 &(HaloTable[currentHalo].SurfP),
		 &(HaloTable[currentHalo].Phi0),
		 &(HaloTable[currentHalo].cNFW) 
		 );
	  IDmap[currentHalo] = HaloTable[currentHalo].ID;
	  HaloTable[currentHalo].AHFID = HaloTable[currentHalo].ID; 
	  HaloTable[currentHalo].SnapID = iFile;
	  HaloTable[currentHalo].nAvatars = 1;
	  HaloTable[currentHalo].AvatarList = calloc(HaloTable[currentHalo].nAvatars,sizeof(MyIDtype));
	  HaloTable[currentHalo].AvatarList[0] = currentHalo;
	  Avatar[currentHalo] = currentHalo;
	  HaloTable[currentHalo].ProgAvatarFlag = 0;
	  HaloTable[currentHalo].TroubleFlag = 0;

	  //printf("%llu : %llu\n", currentHalo,IDmap[currentHalo]);
	  currentHalo++;
	}
      fclose(fp);
      
    }

  if(TotNhalos != (currentHalo))
    {
      fprintf(stderr,"Error: Problem in halo catalogues: %llu,%llu\nExiting...\n",TotNhalos,currentHalo);
      exit(0);
    }
  SubTree = calloc(TotNhalos,sizeof(MyIDtype));
  for(iHalo=0;iHalo<TotNhalos;iHalo++)
    {
      SubTree[iHalo] = IDsearch(HaloTable[iHalo].hostHalo);
      //printf("halo %llu:%llu hosted in %llu: %llu\n", iHalo,IDmap[iHalo], SubTree[iHalo],HaloTable[iHalo].hostHalo);
    }
  
  for(iHalo=0;iHalo<TotNhalos;iHalo++)
    {
      hosthalo = SubTree[iHalo];
      if(hosthalo < NULLPOINT)
	{
	  if(HaloTable[hosthalo].Mvir < 1.e-30)
	    HaloTable[iHalo].Mvir = 0.;
	}
    }


  read_particles(snapnum);
  sprintf(filename,"%s/snapdir_%03d/%s%03d",gadgetfolder,(int)snapnum,gadgetPrefix,(int)snapnum);
  i = (MyIDtype) gadget_load_snapshot(filename,16,P,PIDmap );
  printf("Total particle : %llu\n",i);
}


int gadget_load_snapshot(char *fname, int files, struct Gadget_particle *P, int *PIDmap)
{
  FILE *fd;
  char buf[200];
  int i, j, k, dummy, ntot_withmasses,NumPart,Ngas;
  int t, n, off, pc, pc_new, pc_sph;
  int *Id;
  struct gadget_io_header header1;
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i = 0, pc = 1; i < files; i++, pc = pc_new)
    {
      if(files > 1)
	sprintf(buf, "%s.%d", fname, i);
      else
	sprintf(buf, "%s", fname);

      if(!(fd = fopen(buf, "r")))
	{
	  printf("can't open file `%s`\n", buf);
	  exit(0);
	}

      printf("reading `%s' ...\n", buf);
      fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files == 1)
	{
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    {
	      NumPart += header1.npart[k];
	    }
	  Ngas = header1.npart[0];
	}
      else
	{
	  
	  for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
	    {
	      NumPart += header1.npartTotal[k];
	      printf("Total num %d : %d\n",k,header1.npartTotal[k]);
	    }
	  Ngas = header1.npartTotal[0];
	}

      for(k = 0, ntot_withmasses = 0; k < 6; k++)
	{
	  if(header1.mass[k] == 0)
	    ntot_withmasses += header1.npart[k];
	}

      if(i == 0)
	{
	  printf("Total numpart : %d\n",NumPart);
	  
	  //allocate_memory();
	  P = malloc(NumPart * sizeof(struct Gadget_particle));
	  Id = malloc(NumPart * sizeof(int));
	  PIDmap = malloc(NumPart * sizeof(int));

	  P--;
	  Id--;
	  PIDmap--;
	}
      for(k=0;k<6;k++)
	{
	  printf("N[%d] : %d\n",k,header1.npart[k]);
	}
      for(k=0;k<6;k++)
	{
	  printf("Total N[%d] : %d\n",k,header1.npartTotal[k]);
	}
      for(k=0;k<6;k++)
	{
	  printf("mass[%d] : %lf\n",k,(double)header1.mass[k]);
	}
      printf("time : %lf\n",header1.time);
      printf("redshift : %lf\n",header1.redshift);
      printf("flag_sfr : %d\n", header1.flag_sfr);
      printf("flag_feedback : %d\n", header1.flag_feedback);
      printf("flag_cooling : %d\n", header1.flag_cooling);
      printf("numfiles : %d\n", header1.num_files);
      printf("BoxSize : %lf\n",header1.BoxSize);
      printf("Omega0 : %lf\n",header1.Omega0);
      printf("OmegaLambda : %lf\n",header1.OmegaLambda);
      printf("HubbleParam : %lf\n",header1.HubbleParam);

      for(k=0;k<22;k++)
	{
	  printf("fillers[%d] : %d\n",k, header1.filler[k]);
	}
      SKIP;
      printf("dummy = %d\n",dummy);
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
	      printf("%d => %f\n",pc_new,P[pc_new].Pos[0]);
	      printf("%d => %f\n",pc_new,P[pc_new].Pos[1]);
	      printf("%d => %f\n",pc_new,P[pc_new].Pos[2]);
	      exit(0);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;


      SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      fread(&Id[pc_new], sizeof(int), 1, fd);
	      pc_new++;
	    }
	}
      SKIP;


      if(ntot_withmasses > 0)
	SKIP;
      for(k = 0, pc_new = pc; k < 6; k++)
	{
	  for(n = 0; n < header1.npart[k]; n++)
	    {
	      P[pc_new].Type = k;

	      if(header1.mass[k] == 0)
		fread(&P[pc_new].Mass, sizeof(float), 1, fd);
	      else
		P[pc_new].Mass = header1.mass[k];
	      pc_new++;
	    }
	}
      if(ntot_withmasses > 0)
	SKIP;


      if(header1.npart[0] > 0)
	{
	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
	      fread(&P[pc_sph].U, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  SKIP;
	  for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	    {
	      fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
	      pc_sph++;
	    }
	  SKIP;

	  if(header1.flag_cooling)
	    {
	      SKIP;
	      for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
		{
		  fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
		  pc_sph++;
		}
	      SKIP;
	    }
	  else
	    for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
	      {
		P[pc_sph].Ne = 1.0;
		pc_sph++;
	      }
	}

      fclose(fd);
    }
  printf("%d %f %f %f\n",Id[NumPart],P[NumPart].Pos[0],P[NumPart].Vel[0],P[NumPart].Mass);
  for(i=1;i<=header1.npartTotal[1];i++)
    {
      printf("%d %f %f %f\n",(int)Id[i],P[i].Pos[0],P[i].Vel[0],P[i].Mass);
      //printf("%d => %d\n",i,(int)Id[i]);
      //PIDmap[Id[i]] = i;
    }
  return NumPart;
}

void hbtmaphalos(char file[MAXSTRING])
{
  MyIDtype ihalo,jhalo,j,counthalo,ahf_haloid,countuseable;
  const unsigned long long MAXHBT = 2000000;
  int *original_used;
  char line[MAXSTRING],buff[MAXSTRING];
  FILE* fp;
  
  sprintf(buff,"%s/HBT_added_halos.all",file);
  HBT = calloc(MAXHBT,sizeof(struct HBT_halos));
  fp = fopen(buff, "r");
  fgets(line,MAXSTRING,fp);
  counthalo = 0;
  while((fgets(line,MAXSTRING,fp)) != NULL)
    {
      //printf("%s",line);
      sscanf(line,"%llu %llu %lld %llu %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g",
	     &(HBT[counthalo].SubHaloID),
	     &(HBT[counthalo].HostID),
	     &(HBT[counthalo].AHFID),
	     &(HBT[counthalo].Nbound),
	     &(HBT[counthalo].RvirEstimate),
	     &(HBT[counthalo].Vmax),
	     &(HBT[counthalo].Rmax),
	     &(HBT[counthalo].X),
	     &(HBT[counthalo].Y),
	     &(HBT[counthalo].Z),
	     &(HBT[counthalo].Vx),
	     &(HBT[counthalo].Vy),
	     &(HBT[counthalo].Vz),
	     &(HBT[counthalo].Mvir),
	     &(HBT[counthalo].Rvir),
	     &(HBT[counthalo].Rhalf),
	     &(HBT[counthalo].RPoisson),
	     &(HBT[counthalo].R2Sig),
	     &(HBT[counthalo].Req),
	     &(HBT[counthalo].Rtidal)	     
	     );
      counthalo++;
    }
  
  fclose(fp);
  HBT = (struct HBT_halos *) realloc (HBT, counthalo * sizeof(struct HBT_halos));
  countuseable = 0;

  original_used = calloc(TotNhalos,sizeof(int));
  HBThaloTable = calloc(counthalo,sizeof(struct HALOPROPS));
  for(ihalo=0;ihalo<counthalo;ihalo++)
    {
      if(HBT[ihalo].AHFID >= 0)
	{
	  ahf_haloid = (MyIDtype) HBT[ihalo].AHFID;
	  //force to use HBT data
	  ahf_haloid = NULLPOINT;
	  //ahf_haloid = IDsearch(ahf_haloid);
	  if(ahf_haloid < NULLPOINT)
	    {
	      if(original_used[ahf_haloid] == 0)
		{
		  original_used[ahf_haloid] = 1;
		  HBThaloTable[ihalo].ID = HBT[ihalo].SubHaloID;
		  HBThaloTable[ihalo].SnapID = HaloTable[ahf_haloid].SnapID;
		  HBThaloTable[ihalo].hostHalo = HBT[ihalo].HostID;
		  HBThaloTable[ihalo].Mvir = HaloTable[ahf_haloid].Mvir;
		  HBThaloTable[ihalo].Xc = HaloTable[ahf_haloid].Xc;
		  HBThaloTable[ihalo].Yc = HaloTable[ahf_haloid].Yc;
		  HBThaloTable[ihalo].Zc = HaloTable[ahf_haloid].Zc;
		  HBThaloTable[ihalo].VXc = HaloTable[ahf_haloid].VXc;
		  HBThaloTable[ihalo].VYc = HaloTable[ahf_haloid].VYc;
		  HBThaloTable[ihalo].VZc = HaloTable[ahf_haloid].VZc;	
		  HBThaloTable[ihalo].Vmax = HaloTable[ahf_haloid].Vmax;
		  HBThaloTable[ihalo].Rvir = HaloTable[ahf_haloid].Rvir;
		}
	      else
		{
		  printf("duplicated map %llu\n",HBT[ihalo].SubHaloID);
		  exit(0);
		}
	    }
	  else
	    {
	      HBThaloTable[ihalo].ID = HBT[ihalo].SubHaloID;
	      HBThaloTable[ihalo].SnapID = HBThaloTable[ihalo].ID/MAXHALOPERSNAP;
	      HBThaloTable[ihalo].hostHalo = HBT[ihalo].HostID;
	      if(HBT[ihalo].AHFID > 0)
		HBThaloTable[ihalo].AHFID = HBT[ihalo].AHFID;
	      else
		HBThaloTable[ihalo].AHFID = NULLPOINT;
	      HBThaloTable[ihalo].Mvir = HBT[ihalo].Mvir*HBT2AHFmass;
	      HBThaloTable[ihalo].npart = HBT[ihalo].Nbound;
	      HBThaloTable[ihalo].Xc = HBT[ihalo].X;
	      HBThaloTable[ihalo].Yc = HBT[ihalo].Y;
	      HBThaloTable[ihalo].Zc = HBT[ihalo].Z;
	      HBThaloTable[ihalo].VXc = HBT[ihalo].Vx;
	      HBThaloTable[ihalo].VYc = HBT[ihalo].Vy;
	      HBThaloTable[ihalo].VZc = HBT[ihalo].Vz;
	      HBThaloTable[ihalo].Vmax = HBT[ihalo].Vmax;
	      HBThaloTable[ihalo].Rvir = HBT[ihalo].Rvir;
	      //printf("cannot find halo ID %llu %lld\n",HBT[ihalo].SubHaloID, HBT[ihalo].AHFID);
	      //exit(0);
	    }
	}
      else
	{
	  HBThaloTable[ihalo].ID = HBT[ihalo].SubHaloID;
	  HBThaloTable[ihalo].SnapID = HBThaloTable[ihalo].ID/MAXHALOPERSNAP;
	  HBThaloTable[ihalo].hostHalo = HBT[ihalo].HostID;

	  if(HBT[ihalo].AHFID > 0)
	    HBThaloTable[ihalo].AHFID = HBT[ihalo].AHFID;
	  else
	    HBThaloTable[ihalo].AHFID = NULLPOINT;

	  HBThaloTable[ihalo].Mvir = HBT[ihalo].Mvir*HBT2AHFmass;
	  HBThaloTable[ihalo].npart = HBT[ihalo].Nbound;
	  HBThaloTable[ihalo].Xc = HBT[ihalo].X;
	  HBThaloTable[ihalo].Yc = HBT[ihalo].Y;
	  HBThaloTable[ihalo].Zc = HBT[ihalo].Z;
	  HBThaloTable[ihalo].VXc = HBT[ihalo].Vx;
	  HBThaloTable[ihalo].VYc = HBT[ihalo].Vy;
	  HBThaloTable[ihalo].VZc = HBT[ihalo].Vz;
	  HBThaloTable[ihalo].Vmax = HBT[ihalo].Vmax;
	  HBThaloTable[ihalo].Rvir = HBT[ihalo].Rvir;
	}
    }
  printf("total halo = %llu\n",counthalo);
  //printf("useable halo = %llu\n",countuseable);
  free(HBT);
  free(HaloTable);
  TotNhalos = counthalo;
  HaloTable = HBThaloTable;
  resetIDmap();
  for(ihalo=0;ihalo<TotNhalos;ihalo++)
    {
      for(j=0;j<HaloTable[ihalo].nSubhalos;j++)
	{
	  jhalo = HaloTable[ihalo].SubhaloList[j];
	  HaloTable[ihalo].Mvir += HaloTable[jhalo].oriMvir;
	}
    }
#ifdef READPARTICLE
#endif
}

void get_snap_stats()
{
  MyIDtype ihalo,count_halo;
  int isnap;
  float Vmax, thishalo_vel, thishalo_v2;
  long double sum_v2;
  count_halo = 0;
  for(isnap=FIRSTSNAP; isnap < NSNAPS; isnap++)
    {
      //printf("Calculating statistics for snapshot %d\n",isnap);
      Vmax = 0.;
      sum_v2 = 0.;
      for(ihalo=0;ihalo < SnapNhalos[isnap]; ihalo++)
	{
	  thishalo_v2 = HaloTable[count_halo].VXc*HaloTable[count_halo].VXc 
	    + HaloTable[count_halo].VYc*HaloTable[count_halo].VYc
	    + HaloTable[count_halo].VZc*HaloTable[count_halo].VZc;
	  sum_v2 += thishalo_v2;
	  thishalo_vel = sqrt(thishalo_v2);
	  Vmax = MAX(Vmax,thishalo_vel);
	  count_halo++;
	}
      snap_stats[isnap].Vmax = (float) Vmax;
      snap_stats[isnap].Vrms = (float) (sqrt(sum_v2)/SnapNhalos[isnap]);
      //printf("snap %d => %f\n",isnap,snap_stats[isnap].Vmax);
    }
}


void deleteHalos_v1(char file[MAXSTRING])
{
  FILE* fp;
  MyIDtype selecthaloID,selecthalo,count;
  char line[MAXSTRING];
  printf("Start deleting halos from %s\n",file);
  fp = fopen(file, "r");
  count = 0;
  while((fgets(line,MAXSTRING,fp)) != NULL)
    {
 
      sscanf(line,"%llu",&selecthaloID);
      //printf("read %llu\n",selecthaloID);
      selecthalo = IDsearch(selecthaloID);

      //printf("%llu\t%llu\n",selecthaloID,selecthalo);

      if(selecthalo < NULLPOINT)
	{
	  //printf("%llu\n",selecthaloID);
	  //IDmap[selecthalo] = NULLPOINT;
	  HaloTable[selecthalo].ID = NULLPOINT;
	  
	  count++;
	}
    }
  fclose(fp);
  TotNhalosUsed -= count;
  resetIDmap();
  printf("TotNhalos: %llu\nDeleted Halos: %llu\nTotNhalosUsed: %llu\n\n",TotNhalos,count,TotNhalosUsed);
}
void resetIDmap()
{
  MyIDtype ihalo,jhalo,khalo,hostid,is_exist,counthalo,j,k,l,count,dummyid;
  int isnap,passcheck;
  printf("Start Reset IDmap\n");
  qsort(HaloTable,TotNhalos, sizeof(struct HALOPROPS), compareID);
  count = 0;
  for(isnap = FIRSTSNAP; isnap < NSNAPS; isnap++)
    {
      SnapNhalos[isnap] = 0;
    }
  for(ihalo=0;ihalo<TotNhalos;ihalo++)
    {
      if(HaloTable[ihalo].ID < NULLPOINT)
	{
	  SnapNhalos[HaloTable[ihalo].SnapID]++;
	  free(HaloTable[ihalo].AvatarList);
	  count++;
	}
      else
	break;
    }
  TotNhalos = count;
  TotNavatars = TotNhalos;
  printf("TotNhalos = %llu\n",TotNhalos);
  HaloTable = realloc(HaloTable,TotNhalos*sizeof(struct HALOPROPS));
  free(IDmap);
  free(SubTree);
  free(Avatar);
  IDmap = calloc(TotNhalos,sizeof(MyIDtype));
  SubTree = calloc(TotNhalos,sizeof(MyIDtype));
  Avatar = calloc(TotNhalos,sizeof(MyIDtype));
  for(ihalo=0;ihalo<TotNhalos;ihalo++)
    {
      IDmap[ihalo] = HaloTable[ihalo].ID;
      HaloTable[ihalo].oriMvir = HaloTable[ihalo].Mvir;
      //SubTree[ihalo] = IDsearch(HaloTable[ihalo].hostHalo);
      HaloTable[ihalo].SubhaloList = malloc(0);
      HaloTable[ihalo].nSubhalos = 0;
      HaloTable[ihalo].nAvatars = 1;
      HaloTable[ihalo].AvatarList = calloc(HaloTable[ihalo].nAvatars,sizeof(MyIDtype));
      HaloTable[ihalo].AvatarList[0] = ihalo;
      Avatar[ihalo] = ihalo;
    }
  printf("Finish setup initials\n");
  // a very expensive way to make an inclusive catalogue
  for(ihalo=0;ihalo<TotNhalos;ihalo++)
    {
      SubTree[ihalo] = IDsearch(HaloTable[ihalo].hostHalo);
      hostid = SubTree[ihalo];
      //printf("%llu => %llu\n",ihalo,SubTree[ihalo]);
      if(SubTree[ihalo] < NULLPOINT)
	{
	  HaloTable[hostid].nSubhalos++;
	  HaloTable[hostid].SubhaloList = realloc(HaloTable[hostid].SubhaloList,HaloTable[hostid].nSubhalos*sizeof(MyIDtype));
	  HaloTable[hostid].SubhaloList[HaloTable[hostid].nSubhalos-1] = ihalo;
	}
    }
  printf("Finish setup initials for inclusive mass\n");
  for(ihalo=0;ihalo<TotNhalos;ihalo++)
    {
      //printf("checking halo %llu\n",ihalo);
      passcheck = 0;
      //printf("ihalo = %llu -: %llu subhalos\n", ihalo,HaloTable[ihalo].nSubhalos);
      while(passcheck != 1)
	{
	  passcheck = 1;
	  //printf("round up :\n");
	  for(j=0; j<HaloTable[ihalo].nSubhalos; j++)
	    {
	      //printf("j = %llu :",j);
	      jhalo = HaloTable[ihalo].SubhaloList[j];
	      //printf("%llu: ",jhalo);
	      for(k=0;k<HaloTable[jhalo].nSubhalos;k++)
	  	{
	  	  khalo = HaloTable[jhalo].SubhaloList[k];
	  	  //printf("%llu ",khalo);
	  	  is_exist = 0;
	  	  for(l=0;l<HaloTable[ihalo].nSubhalos;l++)
	  	    {
	  	      if(khalo == HaloTable[ihalo].SubhaloList[l])
	  		{
	  		  is_exist = 1;
	  		}
	  	    }
	  	  if(is_exist == 0)
	  	    {
	  	      passcheck = 0;
	  	      HaloTable[ihalo].nSubhalos++;
	  	      HaloTable[ihalo].SubhaloList = realloc(HaloTable[ihalo].SubhaloList,HaloTable[ihalo].nSubhalos*sizeof(HaloTable[ihalo].nSubhalos));
	  	      HaloTable[ihalo].SubhaloList[HaloTable[ihalo].nSubhalos-1] = khalo;
	  	    }

	  	}

	    }
	  //printf("\n");
	}
    }
  printf("Finish Reset IDmap\n");
}

void addHalo_v1(char file[MAXSTRING])
{
  FILE* fp;
  MyIDtype ihalo,j,k,count,dummyid;
  int dummyint;
  float dummyfloat;
  char line[MAXSTRING];
  printf("Start adding halos from %s\n",file);
  fp = fopen(file,"r");
  count = 0;
  while((fgets(line,MAXSTRING,fp)) != NULL)
    {
      count++;
    }
  fclose(fp);
  ihalo = TotNhalos;
  if(count <= MAXUSEABLE && count > 0)
    {
      TotNhalos += count;
      TotNavatars = TotNhalos;
      TotNhalosUsed = TotNhalos;
      printf("TotNhalos: %llu\nAdded halos: %llu\nTotNhalosUsed: %llu\n\n",TotNhalos,count,TotNhalosUsed);

      HaloTable = realloc(HaloTable, TotNhalos*sizeof(struct HALOPROPS));
      printf("finish realloc hALOTABLE\n");
      Avatar= realloc(Avatar,TotNhalos*sizeof(MyIDtype));
      printf("finish realloc aVATAR\n");
      SubTree = realloc(SubTree,TotNhalos*sizeof(MyIDtype));
      printf("finish realloc SUBTREE\n");
      IDmap = realloc(IDmap,TotNhalos*sizeof(MyIDtype));
      printf("finish realloc IDmap\n");
    }

  fp = fopen(file,"r");
  while((fscanf(fp,"%d %llu %g %g %g %g %g %g %g %g %g %g",
		&(HaloTable[ihalo].SnapID),
		&(HaloTable[ihalo].ID),
		&(HaloTable[ihalo].Mvir),
		&(HaloTable[ihalo].Rvir),
		&(HaloTable[ihalo].r2), 
		&(HaloTable[ihalo].Vmax),	
		&(HaloTable[ihalo].Xc),
		&(HaloTable[ihalo].Yc),
		&(HaloTable[ihalo].Zc),
		&(HaloTable[ihalo].VXc),
		&(HaloTable[ihalo].VYc),
		&(HaloTable[ihalo].VZc)	   
		)) != EOF)
    {
      //printf("%llu %llu\n",ihalo,HaloTable[ihalo].ID);
      IDmap[ihalo] = HaloTable[ihalo].ID;
      HaloTable[ihalo].npart = 20;
      HaloTable[ihalo].nAvatars = 1;
      HaloTable[ihalo].AvatarList = calloc(HaloTable[ihalo].nAvatars,sizeof(MyIDtype));
      HaloTable[ihalo].AvatarList[0] = ihalo;
      HaloTable[ihalo].AHFID = NULLPOINT;
      Avatar[ihalo] = ihalo;
      HaloTable[ihalo].ProgAvatarFlag = 0;
      HaloTable[ihalo].TroubleFlag = 0;
      SubTree[ihalo] = NULLPOINT;
      HaloTable[ihalo].Xc *= Mpc2kpc;
      HaloTable[ihalo].Yc *= Mpc2kpc;
      HaloTable[ihalo].Zc *= Mpc2kpc;
      ihalo++;
    }
  fclose(fp);

  resetIDmap();
  printf("Finish adding halos\n");
}

void load_particles(MyIDtype load_id)
{
  FILE *f;
  char filename[MAXSTRING];
  char dummystr[MAXSTRING];
  float dummyfloat;
  MyIDtype i, j, dummyID,haloid,u_haloid,nhalos,nparts;
  unsigned int snapid;
  struct particle_data dummyhalo;

  

  snapid = (unsigned int) (load_id/MAXHALOPERSNAP);
  strncpy(filename,"",sizeof(filename));
  (void) getFilename(filename,snapid);
  strcat(filename,".AHF_particles");
  f = fopen(filename,"r");
  printf("%s\n",filename);
  if(f == NULL)
    {
      printf("could not open %s\nABORTING\n",filename);
      exit(1);
    }
  //fgets(dummystr,MAXSTRING,f);
  fscanf(f, "%llu", &nhalos);
  if(nhalos != SnapNhalos[snapid])
    {
      printf("Mismatch total halos in snapshot %d\n",snapid);
      exit(0);
    }
  for(i=0; i<SnapNhalos[snapid]; i++)
    {
      fscanf(f, "%llu %llu",&nparts,&u_haloid);
      if(u_haloid == load_id)
	{
	  haloid = IDsearch(u_haloid);
	  //printf("halo %llu => %llu : npart:%llu\n",u_haloid,haloid,nparts);
	  if(nparts != HaloTable[haloid].npart)
	    {
	      printf("Mismatch total particles in halo %llu => %llu\n",u_haloid,haloid);
	      printf("HaloTable :%llu  read:%llu\n",HaloTable[haloid].npart,nparts);
	      exit(0);
	    }
	  HaloTable[haloid].Particles     = (struct particle_data *) calloc(nparts, sizeof(struct particle_data));

	  for(j=0; j<nparts; j++)
	    {
	      //printf("read particle %llu\n",j);
	      fscanf(f, "%llu %g %g %g %g %g %g %g",

		     &(HaloTable[haloid].Particles[j].ParticleID),
		     &(HaloTable[haloid].Particles[j].ParticleEnergy),       // [Msun/h (km/sec)^2]

		     &(HaloTable[haloid].Particles[j].X),                    // [kpc/h]
		     &(HaloTable[haloid].Particles[j].Y),                    // [kpc/h]
		     &(HaloTable[haloid].Particles[j].Z),                    // [kpc/h]
		     &(HaloTable[haloid].Particles[j].Vx),                   // [km/sec]
		     &(HaloTable[haloid].Particles[j].Vy),                   // [km/sec]
		     &(HaloTable[haloid].Particles[j].Vz)                    // [km/sec]		
		     );	  
	    }
	  qsort(HaloTable[haloid].Particles,nparts, sizeof(struct particle_data), compareParticleEnergy);
	}
      else
	{
	  for(j=0; j<nparts; j++)
	    {
	      //printf("read particle %llu\n",j);
	      fscanf(f, "%llu %g %g %g %g %g %g %g",

		     &(dummyhalo.ParticleID),
		     &(dummyhalo.ParticleEnergy),       // [Msun/h (km/sec)^2]

		     &(dummyhalo.X),                    // [kpc/h]
		     &(dummyhalo.Y),                    // [kpc/h]
		     &(dummyhalo.Z),                    // [kpc/h]
		     &(dummyhalo.Vx),                   // [km/sec]
		     &(dummyhalo.Vy),                   // [km/sec]
		     &(dummyhalo.Vz)                    // [km/sec]		
		     );	  
	    }
	}
    }
  fclose(f);
  
}


void read_particles(unsigned int slotid)
{
  FILE *f;
  char filename[MAXSTRING];
  char dummystr[MAXSTRING];
  float dummyfloat;
  MyIDtype i, j, dummyID,haloid,u_haloid,nhalos,nparts;
  unsigned int snapid;
  struct particle_data dummyhalo;

  snapid =  slotid;
  strncpy(filename,"",sizeof(filename));

  if(output.outputFormat > 1.129 && output.outputFormat < 1.131)
    {
      (void) getHBTFilename(filename,snapid);
      strcat(filename,".HBT_particles");
    }
  else
    {
      (void) getFilename(filename,snapid);
      strcat(filename,".AHF_particles");
    }
  f = fopen(filename,"r");
  printf("%s\n",filename);
  if(f == NULL)
    {
      printf("could not open %s\nABORTING\n",filename);
      exit(1);
    }
  //fgets(dummystr,MAXSTRING,f);
  fscanf(f, "%llu", &nhalos);
  if(nhalos != SnapNhalos[snapid])
    {
      printf("Mismatch total halos in snapshot %d\n",snapid);
      exit(0);
    }
  for(i=0; i<SnapNhalos[snapid]; i++)
    {
      fscanf(f, "%llu %llu",&nparts,&u_haloid);
      haloid = IDsearch(u_haloid);
      //printf("halo %llu => %llu : npart:%llu\n",u_haloid,haloid,nparts);
      if(nparts != HaloTable[haloid].npart)
	{
	  printf("Mismatch total particles in halo %llu => %llu\n",u_haloid,haloid);
	  printf("HaloTable :%llu  read:%llu\n",HaloTable[haloid].npart,nparts);
	  exit(0);
	}
      HaloTable[haloid].Particles     = (struct particle_data *) calloc(nparts, sizeof(struct particle_data));

      for(j=0; j<nparts; j++)
	{
	  //printf("read particle %llu\n",j);
	  fscanf(f, "%llu %g %g %g %g %g %g %g",

		 &(HaloTable[haloid].Particles[j].ParticleID),
		 &(HaloTable[haloid].Particles[j].ParticleEnergy),       // [Msun/h (km/sec)^2]

		 &(HaloTable[haloid].Particles[j].X),                    // [kpc/h]
		 &(HaloTable[haloid].Particles[j].Y),                    // [kpc/h]
		 &(HaloTable[haloid].Particles[j].Z),                    // [kpc/h]
		 &(HaloTable[haloid].Particles[j].Vx),                   // [km/sec]
		 &(HaloTable[haloid].Particles[j].Vy),                   // [km/sec]
		 &(HaloTable[haloid].Particles[j].Vz)                    // [km/sec]		
	   );	  
	}
      qsort(HaloTable[haloid].Particles,nparts, sizeof(struct particle_data), compareParticleEnergy);
    }
  fclose(f);
  
}
void read_particles_binary()
{
  FILE *fp_npart,*fp_particle, *fp_uid;
  char filename[MAXSTRING];
  MyIDtype i,j,dummy;
  if(output.outputFormat > 1.129 && output.outputFormat < 1.131)
    {
      sprintf(filename,"%s/%s_particles.bin",FolderName,inputFile);
    }
  else
    {
      sprintf(filename,"%s/all_particles.bin",FolderName);
    }
  fp_particle = fopen(filename, "rb");


  fread (&dummy , 1 , sizeof(MyIDtype) , fp_particle);
  if(dummy != TotNhalos)
    {
      printf("%llu != %llu\n",dummy,TotNhalos);
      exit(0);
    }
  for(i=0;i<TotNhalos;i++)
    {
      //fread (&(HaloTable[i].ID) , 1 , sizeof(MyIDtype) , fp_uid );
      //fread (&(HaloTable[i].npart) , 1 , sizeof(MyIDtype) , fp_npart );
      HaloTable[i].Particles = calloc(HaloTable[i].npart,sizeof(struct particle_data));
      for(j=0;j<HaloTable[i].npart;j++)
	{
	  fread (&(HaloTable[i].Particles[j].ParticleID) , 1 , sizeof(MyIDtype), fp_particle);
	}
    }
  fclose(fp_particle);
}

void getFilename(char* filename,unsigned int snapnum)
{
  char snapstr[MAXSTRING];
  int len;
  struct dirent *pDirent;
  DIR *pDir;
  unsigned int iFile, i;
  MyIDtype currentHalo, nHalo;
  char keyword[MAXSTRING];
  char dummystr[MAXSTRING];
  char zstr[MAXSTRING];
  char *returnstr,*finalstr; 
  //printf("getFilename %d\n",snapnum);
  sprintf(keyword,"%s%03d.",FilePrefix,snapnum);
  pDir = opendir(FolderName);
  if (pDir == NULL) 
    {
      printf ("Cannot open directory '%s'\n", FolderName);
      exit(0);
    }
  while ((pDirent = readdir(pDir)) != NULL) {
    if((returnstr=strstr(pDirent->d_name,keyword)))
      {
	len = strlen(returnstr);
	sprintf(dummystr,"%s",returnstr+strlen(keyword));
	len = len - strlen(keyword);
	returnstr = strstr(dummystr,".AHF_");
	len = len - strlen(returnstr);
	strncpy(zstr,"",sizeof(zstr));
	strncpy(zstr,dummystr,len);
	//printf("z=%s\n",zstr);
	sprintf(filename,"%s/%s%s",FolderName,keyword,zstr);
	break;
      }
  } 
  closedir (pDir);
}

void getHBTFilename(char* filename,unsigned int snapnum)
{
  char snapstr[MAXSTRING];
  int len;
  struct dirent *pDirent;
  DIR *pDir;
  unsigned int iFile, i;
  MyIDtype currentHalo, nHalo;
  char keyword[MAXSTRING];
  char dummystr[MAXSTRING];
  char zstr[MAXSTRING];
  char *returnstr,*finalstr; 
  //printf("getFilename %d\n",snapnum);
  sprintf(keyword,"%s%03d.",HBTFilePrefix,snapnum);
  pDir = opendir(HBTFolderName);
  if (pDir == NULL) 
    {
      printf ("Cannot open directory '%s'\n", HBTFolderName);
      exit(0);
    }
  while ((pDirent = readdir(pDir)) != NULL) {
    if((returnstr=strstr(pDirent->d_name,keyword)))
      {
	len = strlen(returnstr);
	sprintf(dummystr,"%s",returnstr+strlen(keyword));
	len = len - strlen(keyword);
	returnstr = strstr(dummystr,".HBT_");
	len = len - strlen(returnstr);
	strncpy(zstr,"",sizeof(zstr));
	strncpy(zstr,dummystr,len);
	//printf("z=%s\n",zstr);
	sprintf(filename,"%s/%s%s",HBTFolderName,keyword,zstr);
	break;
      }
  } 
  closedir (pDir);
}

