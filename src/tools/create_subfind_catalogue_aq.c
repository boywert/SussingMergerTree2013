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


#define MINCOMMON      10   // we only cross-correlate haloes if they at least share MINCOMMON particles
#define NSATinHOSTfrac 0.5  // at least NSATinHOSTfrac the subhaloes particles need to be inside the prospective host

#define SO_VEL_DISPERSIONS

//#define DEBUG
//#define MTREE_SELF

// writes output that readily allows to find mergers (note, the _idx file will be replaced by _merger)
//#define MERGER_RATIO   0.25

#define CLUES_WM3

#define NSNAPS 62
#define MAXSTRING 4096
#define MaxHiResID 18535972

static float default_float = 2e38; 
static float Mpc2kpc = 1000.;
static float kpc2km = 3.08567758e16;
static float GagetUnit2Msun = 1e10;

typedef unsigned int MyIDType;

struct io_header
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
};

struct halo_catalogue
{
  long long TotNids;
  MyIDType *SubOffset;
  MyIDType *IdList;
  float *IdPos;
  float *IdVel;
  float *IdBindingEgy;
  int TotNsubhalos;
  int TotNgroups;
  int *SubLen;
  int *SubParentHalo;
  int *IdToHalo;
  float *Group_M_Crit200;
  float *Group_R_Crit200;
  int *GroupNsubs;
  int *GroupLen;
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
} CatA,CatB;

struct ahf_halos 
{
  unsigned long long  ID;
  unsigned long long   hostHalo;
  MyIDType  numSubStruct;
  float Mvir;
  MyIDType npart;
  float Xc, Yc, Zc;
  float VXc, VYc, VZc;
  float Rvir, Rmax;
  float r2;
  float mbp_offset, com_offset;
  float Vmax, v_esc, sigV;
  float lambda, lambdaE;
  float Lx, Ly, Lz;
  float b,c;
  float Eax, Eay, Eaz;
  float Ebx, Eby, Ebz;
  float Ecx, Ecy, Ecz;
  float ovdens;
  float nbins;
  float fMhires;
  float Ekin, Epot;
  float SurfP;
  float Phi0;
  float cNFW; 
} ahf_id; 


long    *main_progeniter;
long    nHalos[2];
long    PidMax=-1;
long    PidMin=1234567890;
double *merging_rate, *delta_out, *forming_rate;

long unsigned  nlines;       // number of lines in *_mtree file
long unsigned  nstree;       // for how many halos to write *_stree file
long          *ihost;        // assigns one host to each halo

  
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void load_subhalo_catalogue(int num, struct halo_catalogue *cat);
char  OutputDir[MAXSTRING];
void *mymalloc(size_t n);
void myfree(void *ptr);
void long_to_str(char *s, long long n);
void reassign_ids(MyIDType N, MyIDType * ids);

long total_halos(char filename[MAXSTRING]);
int lost_particles(long prog_halo, long merged_halo, int step);
int read_particles(char filename[MAXSTRING], int isimu);
int particle_halo_mapping(int isimu);

float z_list[NSNAPS];
float a[NSNAPS];
char SnapTimeFile[MAXSTRING];
char OUTfolder[MAXSTRING];
void getSnapTime();
double get_z_gadget(int snapid);

int main(int argc, char **argv)
{
  long long i,j,k,l,uniqueID,maxhalopersnap,pid,hosthalo;
  int *subpergroup,numsub,snapid;
  char header_out[4096];
  char h_out[MAXSTRING], p_out[MAXSTRING];
  FILE *fp1,*fp2,*fp3;
  double normalise,factor;
  long long groupNsubs,groupNids,groupfirstsub,groupfirstid;
  double a_c,z;
  strcpy(SnapTimeFile,"data_snaplist.txt");
  strcpy(OUTfolder,"/mnt/lustre/scratch/cs390/SUSSING2013_DATA/datasetII_spin/");
  // sprintf(OUTfolder,"/gpfs/data/Millgas/cs390/SUSSING2013/datasetIII/");
  // use this with datalist_snap
  getSnapTime();

  sscanf(argv[1],"%d",&snapid);
  sprintf(OutputDir,"/mnt/lustre/scratch/cs390/SUSSING2013_DATA/raw_subfind/");
  z = get_z_gadget(snapid);
  a_c = 1./(1+z);
  a[snapid] = (float)a_c;
  z_list[snapid] = (float)z;
  fp3 = fopen("snaplist.txt","a");
  fprintf(fp3,"%d\t%lf\t%lf\n",snapid,a_c,z);
  fclose(fp3);
  //printf("%2.3f\n", z_list[snapid]);
  //sprintf(OutputDir,"/gpfs/data/aquarius/halo_data/Aq-A/4/");
  //snapid = 8;
  maxhalopersnap = pow(10,12);
  load_subhalo_catalogue(snapid, &CatA);


  printf("%d\n", CatA.TotNgroups);

  sprintf(h_out,"%s/datasetII_%03d.z%2.3f.AHF_halos",OUTfolder, snapid, z_list[snapid]);
  sprintf(p_out,"%s/datasetII_%03d.z%2.3f.AHF_particles",OUTfolder, snapid, z_list[snapid]);
  fp1 = fopen(h_out,"w+");
  fp2 = fopen(p_out,"w+");
  k = 0;
  
  subpergroup = calloc(CatA.TotNgroups,sizeof(int));
  sprintf(header_out, 
	  "#ID(1)     hostHalo(2)     numSubStruct(3) Mvir(4) npart(5)        Xc(6)   Yc(7)   Zc(8)   VXc(9)  VYc(10) VZc(11) Rvir(12)        Rmax(13)        r2(14)  mbp_offset(15)  com_offset(16)  Vmax(17)        v_esc(18)       sigV(19)        lambda(20)      lambdaE(21)     Lx(22)  Ly(23)  Lz(24)  b(25)   c(26)   Eax(27) Eay(28) Eaz(29) Ebx(30) Eby(31) Ebz(32) Ecx(33) Ecy(34) Ecz(35) ovdens(36)      nbins(37)       fMhires(38)     Ekin(39)        Epot(40)        SurfP(41)       Phi0(42)        cNFW(43)");
  //printf("%s\n",header_out);
  fprintf(fp1,"%s\n", header_out);
  fprintf(fp2,"%d\n", CatA.TotNsubhalos);
  for(i=0;i<CatA.TotNgroups;i++)
    {
      if(CatA.Group_M_Crit200[i] == 0)
	printf("Mvir :%f Rvir :%f npart :%d\n",CatA.Group_M_Crit200[i],CatA.Group_R_Crit200[i],CatA.GroupLen[i]);
      else
	factor = CatA.Group_M_Crit200[i]/pow(CatA.Group_R_Crit200[i],3);
      for(j=0;j<CatA.GroupNsubs[i]; j++)
	{
	  if(k != CatA.GroupFirstSub[i]+j )
	    exit(0);
	  ahf_id.ID = snapid*maxhalopersnap + k + 1;
	  
	  if(k == CatA.GroupFirstSub[i])
	    ahf_id.hostHalo = 0;
	  else
	    ahf_id.hostHalo = snapid*maxhalopersnap + CatA.GroupFirstSub[i] + 1;

	  if(k == CatA.GroupFirstSub[i])
	    ahf_id.numSubStruct = CatA.GroupNsubs[i];
	  else
	    ahf_id.numSubStruct = 0;

	  if(k == CatA.GroupFirstSub[i])
	    ahf_id.Mvir = CatA.Group_M_Crit200[i] * GagetUnit2Msun;
	  else
	    ahf_id.Mvir = CatA.SubhaloMass[k] * GagetUnit2Msun;
	  //	  if(k == CatA.GroupFirstSub[i])
	  //  printf("halo %ld: M:Mvir = %f Mvir:Rvir^3 = %f\n",k,CatA.Group_M_Mean200[i], pow(CatA.Group_R_Mean200[i],3));
	  ahf_id.npart = CatA.SubLen[k];
	  ahf_id.Xc = CatA.SubhaloPos[3*k]*Mpc2kpc;
	  ahf_id.Yc = CatA.SubhaloPos[3*k+1]*Mpc2kpc;
	  ahf_id.Zc = CatA.SubhaloPos[3*k+2]*Mpc2kpc;
	  ahf_id.VXc = CatA.SubhaloVel[3*k];
	  ahf_id.VYc = CatA.SubhaloVel[3*k+1];
	  ahf_id.VZc = CatA.SubhaloVel[3*k+2];
	  
	  if(k == CatA.GroupFirstSub[i])
	    {
	      printf("");
	      ahf_id.Rvir = CatA.Group_R_Crit200[i]*Mpc2kpc;
	    }
	  else
	    {
	      ahf_id.Rvir = pow(CatA.SubhaloMass[k]/factor,1./3.)*Mpc2kpc;
	    }
	  
	  ahf_id.Rmax = CatA.SubhaloVmaxRad[k]*Mpc2kpc;

       	  ahf_id.r2 = default_float;
	  
	  ahf_id.mbp_offset = default_float;
	  ahf_id.com_offset = default_float;
	  ahf_id.Vmax = CatA.SubhaloVmax[k];
	  ahf_id.v_esc = default_float;
	  ahf_id.sigV = CatA.SubhaloVelDisp[k];
	  ahf_id.lambdaE = default_float;
	  ahf_id.lambda = default_float;

	  normalise = sqrt(ahf_id.Lx*ahf_id.Lx + ahf_id.Ly*ahf_id.Ly + ahf_id.Lz*ahf_id.Lz);

	  double G = 6.67384e-11; // m^3/(kgs^2
	  double m2kpc = 1./3.08567758e19;
	  double m2km = 0.001;
	  double kg2Msun = 1./1.989e30;

	  G *= m2kpc*pow(m2km,2.)/(kg2Msun);
	  ahf_id.lambda = normalise / ahf_id.Mvir / sqrt(2. * G * ahf_id.Mvir * ahf_id.Rvir);

	  /* ahf_id.Lx = default_float; */
	  /* ahf_id.Ly = default_float; */
	  /* ahf_id.Lz = default_float; */
	  
	  
	  ahf_id.Lx /= normalise;
	  ahf_id.Ly /= normalise;
	  ahf_id.Lz /= normalise;

	  
	  
	  ahf_id.b = default_float;
	  ahf_id.c = default_float;
	  ahf_id.Eax = default_float;
	  ahf_id.Eay = default_float;
	  ahf_id.Eaz = default_float;
	  ahf_id.Ebx = default_float;
	  ahf_id.Eby = default_float;
	  ahf_id.Ebz = default_float;
	  ahf_id.Ecx = default_float;
	  ahf_id.Ecy = default_float;
	  ahf_id.Ecz = default_float;
	  ahf_id.ovdens = default_float;
	  ahf_id.nbins = default_float;
	  ahf_id.fMhires = default_float;
	  ahf_id.Ekin = default_float;
	  ahf_id.Epot = default_float;
	  ahf_id.SurfP = default_float;
	  ahf_id.Phi0 = default_float;
	  ahf_id.cNFW = default_float; 

	  fprintf(fp1,"%llu\t%llu\t%lu\t%.8g\t%lu\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\n",
		  ahf_id.ID,
		  ahf_id.hostHalo,
		  ahf_id.numSubStruct,
		  ahf_id.Mvir,
		  ahf_id.npart,
		  ahf_id.Xc,
		  ahf_id.Yc,
		  ahf_id.Zc,
		  ahf_id.VXc,
		  ahf_id.VYc,
		  ahf_id.VZc,
		  ahf_id.Rvir,
		  ahf_id.Rmax,
		  ahf_id.r2,
		  ahf_id.mbp_offset,
		  ahf_id.com_offset,
		  ahf_id.Vmax,
		  ahf_id.v_esc,
		  ahf_id.sigV,
		  ahf_id.lambda,
		  ahf_id.lambdaE,
		  ahf_id.Lx,
		  ahf_id.Ly,
		  ahf_id.Lz,
		  ahf_id.b,
		  ahf_id.c,
		  ahf_id.Eax,
		  ahf_id.Eay,
		  ahf_id.Eaz,
		  ahf_id.Ebx,
		  ahf_id.Eby,
		  ahf_id.Ebz,
		  ahf_id.Ecx,
		  ahf_id.Ecy,
		  ahf_id.Ecz,
		  ahf_id.ovdens,
		  ahf_id.nbins,
		  ahf_id.fMhires,
		  ahf_id.Ekin,
		  ahf_id.Epot,
		  ahf_id.SurfP,
		  ahf_id.Phi0,
		  ahf_id.cNFW 
		  );

	  fprintf(fp2,"%d\t%llu\n",CatA.SubLen[k],ahf_id.ID);	  
	  for(l=0;l<CatA.SubLen[k];l++)
	    {
	      pid = CatA.SubOffset[k]+l;
#ifndef SUBFIND_SAVE_PARTICLELISTS
	      fprintf(fp2,"%lu\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\n",
		      CatA.IdList[pid],
		      default_float,
		      default_float,
		      default_float,
		      default_float,
		      default_float,
		      default_float,
		      default_float
		      );
#else
	      fprintf(fp2,"%lu\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\n",
		      CatA.IdList[pid],
		      CatA.IdBindingEgy[pid],
		      CatA.IdPos[3*pid]*Mpc2kpc,
		      CatA.IdPos[3*pid+1]*Mpc2kpc,
		      CatA.IdPos[3*pid+2]*Mpc2kpc,
		      CatA.IdVel[3*pid]*sqrt(a[snapid]),
		      CatA.IdVel[3*pid+1]*sqrt(a[snapid]),
		      CatA.IdVel[3*pid+2]*sqrt(a[snapid])
		      );
#endif
	    }
	  k++;
	}
    }
}


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
	     &(a[i]),
	     &(z_list[i]),
	     &dummyfloat,
	     &dummyfloat
	     );
      printf("z[%d] = %f\n",i, z_list[i]);
    }
  
}

double get_z_gadget(int snapid)
{
  char filename[MAXSTRING];
  struct io_header header;
  int dummy;
  FILE* fd;
  sprintf(filename,"%s/snapdir_%03d/snap_C02_400_%03d.0",OutputDir,snapid,snapid);
  if(!(fd=fopen(filename,"r")))  //error check and complain
    {
      printf("gadgetIO module can't open file: '%s'\n",filename);
      fflush(stdout);
      exit(0);
    }

  printf("reading: '%s'\n",filename); 
  fflush(stdout);
      
  fread(&dummy, sizeof(dummy), 1, fd);
  fread(&header, sizeof(header), 1, fd);
  fread(&dummy, sizeof(dummy), 1, fd);
  printf("a = %lf\n", header.redshift);
  //exit(0);
  return   header.redshift;
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
	  cat->IdPos =  mymalloc(3 * sizeof(float) * cat->TotNids);
	  cat->IdVel =  mymalloc(3 * sizeof(float) * cat->TotNids);
	  cat->IdBindingEgy =  mymalloc(sizeof(float) * cat->TotNids);

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
	  cat->GroupLen = mymalloc(sizeof(int) * cat->TotNgroups);
	  cat->Group_M_Crit200 = mymalloc(sizeof(float) * cat->TotNgroups);
	  cat->Group_R_Crit200 = mymalloc(sizeof(float) * cat->TotNgroups);
	  cat->GroupFirstSub = mymalloc(sizeof(int) * cat->TotNgroups);
	  cat->Descendant = mymalloc(sizeof(struct descendant_data) * cat->TotNsubhalos);
	  cat->CountProgenitors = mymalloc(sizeof(int) * cat->TotNsubhalos);
	}

      //fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupLen  */
      my_fread(&cat->GroupLen[groupcount], sizeof(int), ngroups, fd);
      fseek(fd, sizeof(int) * ngroups, SEEK_CUR);	/* skip  GroupOffset  */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  GroupMass  */
      fseek(fd, 3 * sizeof(float) * ngroups, SEEK_CUR);	/* skip  GroupPos */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_M_Mean200 */
      fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_R_Mean200 */
      my_fread(&cat->Group_M_Crit200[groupcount], sizeof(float), ngroups, fd);
      my_fread(&cat->Group_R_Crit200[groupcount], sizeof(float), ngroups, fd);
      // fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_M_Crit200 */
      // fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/* skip  Group_R_Crit200 */
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
      my_fread(&cat->SubhaloSpin[3 * subcount], 3 * sizeof(float), nsubhalos, fd);
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

  idcount = 0;
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

#ifdef  SUBFIND_SAVE_PARTICLELISTS
  idcount = 0;
  for(i = 0; i < nFiles; i++)
    {
      sprintf(buf, "%s/groups_%03d/subhalo_posvel_%03d.%d", OutputDir, num, num, i);
      printf("%s/groups_%03d/subhalo_posvel_%03d.%d\n", OutputDir, num, num, i);
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
 
      fseek(fd, sizeof(double), SEEK_CUR);
      
      printf("read all\n");

      my_fread(&cat->IdPos[3 * idcount], 3 * sizeof(float), nids, fd);
      my_fread(&cat->IdVel[3 * idcount], 3 * sizeof(float), nids, fd);
      fseek(fd, sizeof(int)*nids, SEEK_CUR);
      my_fread(&cat->IdBindingEgy[idcount], sizeof(float), nids, fd);
      fclose(fd);

      idcount += nids;
    }
#endif //SUBFIND_SAVE_PARTICLELISTS

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
