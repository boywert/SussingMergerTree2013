#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

int MAXSTRING=1024;

typedef struct HALOPROPS *HALOPROPptr;
typedef struct HALOPROPS
{
  long unsigned  ID;
  long   hostHalo;
  long unsigned  numSubStruct;
  float Mvir;
  long unsigned npart;
  float         Xc, Yc, Zc;
  float         VXc, VYc, VZc;
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
} HALOPROPS;

HALOPROPptr haloprofile[62];
long unsigned TotNhalos[62];

int main(void)
{
  long i;
  int snapshot;
  char file1[MAXSTRING];
  char file2[MAXSTRING];
  
  sprintf(file1,"/export/virgo/Boyd/Millgas/62.5_ALEX/62.6_dm_061.z0.000.AHF_particles");
  sprintf(file2,"/export/virgo/Boyd/Millgas/62.5_ALEX/62.6_dm_061.z0.000.AHF_halos");
  snapshot = 0;
  read_ahf_profile(file1,file2,snapshot);
  for(i=0;i<TotNhalos[snapshot];i++)
    {
      printf("%ld %ld %ld %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
	      haloprofile[snapshot][i].ID,
	      haloprofile[snapshot][i].hostHalo,
	      haloprofile[snapshot][i].numSubStruct,
	      haloprofile[snapshot][i].Mvir,
	      haloprofile[snapshot][i].npart,
	      haloprofile[snapshot][i].Xc,
	      haloprofile[snapshot][i].Yc,
	      haloprofile[snapshot][i].Zc,
	      haloprofile[snapshot][i].VXc,
	      haloprofile[snapshot][i].VYc,
	      haloprofile[snapshot][i].VZc,
	      haloprofile[snapshot][i].Rvir, 
	      haloprofile[snapshot][i].Rmax,
	      haloprofile[snapshot][i].r2,
	      haloprofile[snapshot][i].mbp_offset,
	      haloprofile[snapshot][i].com_offset,
	      haloprofile[snapshot][i].Vmax,
	      haloprofile[snapshot][i].v_esc,
	      haloprofile[snapshot][i].sigV,
	      haloprofile[snapshot][i].lambda,
	      haloprofile[snapshot][i].lambdaE,
	      haloprofile[snapshot][i].Lx,
	      haloprofile[snapshot][i].Ly,
	      haloprofile[snapshot][i].Lx,
	      haloprofile[snapshot][i].b,
	      haloprofile[snapshot][i].c,
	      haloprofile[snapshot][i].Eax,
	      haloprofile[snapshot][i].Eay,
	      haloprofile[snapshot][i].Eaz,
	      haloprofile[snapshot][i].Ebx, 
	      haloprofile[snapshot][i].Eby,
	      haloprofile[snapshot][i].Ebz,
	      haloprofile[snapshot][i].Ecx, 
	      haloprofile[snapshot][i].Ecy, 
	      haloprofile[snapshot][i].Ecy,
	      haloprofile[snapshot][i].ovdens,
	      haloprofile[snapshot][i].nbins,
	      haloprofile[snapshot][i].fMhires,
	      haloprofile[snapshot][i].Ekin,
	      haloprofile[snapshot][i].Epot,
	      haloprofile[snapshot][i].SurfP,
	      haloprofile[snapshot][i].Phi0,
	      haloprofile[snapshot][i].cNFW 
	      );
    }
  
  return 0;
}

int read_ahf_profile(char particlefile[MAXSTRING],char profilefile[MAXSTRING], int snapshot)
{
  FILE *fpin;
  char line[MAXSTRING];
  long i;
  fprintf(stderr,"Reading file %s ...\n", particlefile);
  
  fpin = fopen(particlefile,"r");
  
  if(fpin == NULL)
    {
      fprintf(stderr,"Could not open file %s\nexiting!\n",particlefile);
      exit(0);
    }

  //read number of particle;
  fgets(line,MAXSTRING,fpin);
  sscanf(line, "%ld", &TotNhalos[snapshot]);
  fclose(fpin);
  haloprofile[snapshot] = (HALOPROPptr) calloc(TotNhalos[snapshot], sizeof(HALOPROPS));
  fpin = fopen(profilefile,"r");
  
  if(fpin == NULL)
    {
      fprintf(stderr,"Could not open file %s\nexiting!\n",profilefile);
      exit(0);
    }
  // read header line
  fgets(line,MAXSTRING,fpin);
  for(i=0;i<TotNhalos[snapshot];i++)
    {
      fscanf(fpin,"%ld %ld %ld %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g",
	     &(haloprofile[snapshot][i].ID),
	     &(haloprofile[snapshot][i].hostHalo),
	     &(haloprofile[snapshot][i].numSubStruct),
	     &(haloprofile[snapshot][i].Mvir),
	     &(haloprofile[snapshot][i].npart),
	     &(haloprofile[snapshot][i].Xc),
	     &(haloprofile[snapshot][i].Yc),
	     &(haloprofile[snapshot][i].Zc),
	     &(haloprofile[snapshot][i].VXc),
	     &(haloprofile[snapshot][i].VYc),
	     &(haloprofile[snapshot][i].VZc),
	     &(haloprofile[snapshot][i].Rvir), 
	     &(haloprofile[snapshot][i].Rmax),
	     &(haloprofile[snapshot][i].r2),
	     &(haloprofile[snapshot][i].mbp_offset),
	     &(haloprofile[snapshot][i].com_offset),
	     &(haloprofile[snapshot][i].Vmax),
	     &(haloprofile[snapshot][i].v_esc),
	     &(haloprofile[snapshot][i].sigV),
	     &(haloprofile[snapshot][i].lambda),
	     &(haloprofile[snapshot][i].lambdaE),
	     &(haloprofile[snapshot][i].Lx),
	     &(haloprofile[snapshot][i].Ly),
	     &(haloprofile[snapshot][i].Lx),
	     &(haloprofile[snapshot][i].b),
	     &(haloprofile[snapshot][i].c),
	     &(haloprofile[snapshot][i].Eax),
	     &(haloprofile[snapshot][i].Eay),
	     &(haloprofile[snapshot][i].Eaz),
	     &(haloprofile[snapshot][i].Ebx), 
	     &(haloprofile[snapshot][i].Eby),
	     &(haloprofile[snapshot][i].Ebz),
	     &(haloprofile[snapshot][i].Ecx), 
	     &(haloprofile[snapshot][i].Ecy), 
	     &(haloprofile[snapshot][i].Ecy),
	     &(haloprofile[snapshot][i].ovdens),
	     &(haloprofile[snapshot][i].nbins),
	     &(haloprofile[snapshot][i].fMhires),
	     &(haloprofile[snapshot][i].Ekin),
	     &(haloprofile[snapshot][i].Epot),
	     &(haloprofile[snapshot][i].SurfP),
	     &(haloprofile[snapshot][i].Phi0),
	     &(haloprofile[snapshot][i].cNFW) 
	     );
    }
  fclose(fpin); 
  return 0;
}

