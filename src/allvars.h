#ifndef ALLVARS_H
#define ALLVARS_H

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
#include <dirent.h>
#include <regex.h>

typedef unsigned long long MyIDtype;


static MyIDtype MAXUSEABLE = 18446744073709551614;
static MyIDtype NULLPOINT =  18446744073709551615;

static MyIDtype MAXHALOPERSNAP = 1000000000000;
static MyIDtype box_particle = 270*270*270;
static double kpc2km = 3.08567758e16;
static double year2second = 365.25*24*3600;

static float GagetUnit2Msun = 1e10;

static double p_mass = 9.36395e+08;
static double h = 0.704;
static double Mpc2kpc = 1000.;
static double HBT2AHFmass = 1.e10;
static float default_float = 0;
static int default_int = 0;
static double Boxsize= 62500.;

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

//configuration


#define MAXSTRING 1024

#ifdef AQUARIUS

#define FIRSTSNAP 0
#define LASTSNAP 1023
#define NSNAPS 1024

#else

#define FIRSTSNAP 0
#define LASTSNAP 61
#define NSNAPS 62

#endif
struct progtable 
{
  MyIDtype haloID;
  MyIDtype nProgs;
  MyIDtype *progID;
};

struct outputRAW
{
  float outputFormat;
  char filename[MAXSTRING];
  struct tm* timestamp;
  char Descriptions[MAXSTRING];
  unsigned int TotNsnapshots;
  unsigned int FirstSnapshot;
  unsigned int LastSnapshot;
  MyIDtype nHalos;
  struct progtable *progs;
};

struct HALOPROPS
{
  MyIDtype  ID;
  MyIDtype AHFID;
  unsigned int SnapID;
  MyIDtype   hostHalo;
  MyIDtype  numSubStruct;
  MyIDtype nSubhalos;
  MyIDtype *SubhaloList;
  MyIDtype  nAvatars;
  MyIDtype *AvatarList;
  int ProgAvatarFlag;
  int TroubleFlag;
  float Mvir;
  float oriMvir;
  MyIDtype npart;
  struct particle_data *Particles;
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
};

struct AvatarTable
{
  MyIDtype MainAvatar;
  MyIDtype MainProg;
  MyIDtype NextProg;
  MyIDtype MainProgMerging;
};

struct SNAPSHOT_STATS 
{
  float Vmax;
  float Vrms;
  float *radbin, *xi;
  float *massbin, *massfunction;
};


struct particle_data
{
  MyIDtype ParticleID;
  float ParticleEnergy;
  float X,Y,Z;
  float Vx,Vy,Vz;
};

struct gadget_io_header
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  int  hashtabsize;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 11];	/* fills to 256 Bytes */
};

struct Gadget_particle
{
  float Pos[3];
  float Vel[3];
  float Mass;
  int Type;

  float Rho, U, Temp, Ne;
};

struct subfind_data 
{
  long long TotNids;
  int *SubOffset;
  unsigned int *IdList;
  float *IdPos;
  float *IdVel;
  float *IdBindingEgy;
  int TotNsubhalos;
  int TotNgroups;
  int *SubLen;
  int *SubParentHalo;
  int *IdToHalo;
  unsigned int *GroupOffset;
  float *GroupMass;
  float *GroupPos;
  float *Group_M_Mean200;
  float *Group_R_Mean200;
  float *Group_M_Crit200;
  float *Group_R_Crit200;
  float *Group_M_TopHat200;
  float *Group_R_TopHat200;
  int *GroupContaminationCount;
  float *GroupContaminationMass;
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
  unsigned int *SubhaloMostBoundID;
  int *SubhaloGrNr;  

};

struct HBT_halos
{
  MyIDtype SubHaloID;
  unsigned int SnapID;
  MyIDtype HostID;
  long long AHFID;
  MyIDtype Nbound;
  float RvirEstimate;
  float Vmax;
  float Rmax;
  float X;
  float Y;
  float Z;
  float Vx;
  float Vy;
  float Vz;
  float Mvir;
  float Rvir;
  float Rhalf;
  float RPoisson;
  float R2Sig;
  float Req;
  float Rtidal;
};
extern struct HBT_halos *HBT;


extern struct SNAPSHOT_STATS snap_stats[NSNAPS];
extern void get_snap_stats();

extern struct outputRAW output;

extern MyIDtype TotNavatars;
extern MyIDtype *Avatar;

extern char delFile[MAXSTRING];
extern char addFile[MAXSTRING];
extern char SnapTimeFile[MAXSTRING];

extern char FilePrefix[MAXSTRING];
extern char HBTFilePrefix[MAXSTRING];
extern char FolderName[MAXSTRING];
extern char HBTFolderName[MAXSTRING];
extern char TreeFolder[MAXSTRING];
extern char ppFolder[MAXSTRING];
extern char gadgetfolder[MAXSTRING];
extern char gadgetPrefix[MAXSTRING];
extern MyIDtype SnapNhalos[NSNAPS];
extern float snapTime[NSNAPS];
extern float snapTimeYear[NSNAPS];
extern float expansion_factor[NSNAPS];
extern MyIDtype TotNhalos;
extern MyIDtype TotNhalosUsed;

extern MyIDtype *IDmap;
extern MyIDtype *SubTree;

extern struct HALOPROPS *HaloTable;
extern char inputFile[MAXSTRING];

extern void calculateMassFunction();
extern void merger_analysis(float minmass, float maxmass, int highlim_npart);
extern void snapshot_stats(int snapid, float minmass, float maxmass, int highlim_npart);
extern void main_branch_analysis(float minmass, float maxmass, int highlim_npart);
extern void count_mergers(MyIDtype haloid, int nStep, double binsize,float minmass, float maxmass);
extern void getSnapTime();
extern void makeIDmap();
extern void deleteHalos_v1(char file[MAXSTRING]);
extern void addHalo_v1(char file[MAXSTRING]);
extern void resetIDmap();
extern void getFilename(char* filename,unsigned int snapnum);
extern void getHBTFilename(char* filename,unsigned int snapnum);
extern void load_particles(MyIDtype load_id);
extern void read_particles(unsigned int slotid);
extern void hbtmaphalos(char file[MAXSTRING]);
extern void read_singlesnap(unsigned int snapnum);
extern int gadget_load_snapshot(char *fname, int files, struct Gadget_particle *P, int *PIDmap);
extern int readRawOutputs(char file[MAXSTRING]);
extern MyIDtype assignAvatarMain();
extern MyIDtype assignAvatar();

 
extern int compareParticleEnergy(const void *v1, const void *v2);
extern int compareHaloID(const void *v1, const void *v2);
extern int compareID(const void *v1, const void *v2);
extern int compareMyIDType(const void *v1, const void *v2);
extern MyIDtype IDsearch( MyIDtype searchKey );
extern MyIDtype Generalsearch( MyIDtype searchKey,MyIDtype n_array ,const void *Array  );

extern MyIDtype findtoplevel(MyIDtype hid);
extern void makesubfindout();
extern void printoutprofile();
extern void printouttrees();
extern void printoutparticles_binary();

#endif
