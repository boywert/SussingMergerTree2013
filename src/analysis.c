#include "allvars.h"

MyIDtype *nMergers, *lostmass, *lostmass_host;
MyIDtype **pp_lostmass;
MyIDtype **pp_nmerger;
void calculateMassFunction()
{
  MyIDtype ihalo, *c_massfunction, usedhalos;
  double *massfunction,*m;
  double Mass,minmass,maxmass,logbin;
  double logminmass, logmaxmass;
  int i,n_bin,block;
  minmass = 1.;
  maxmass = 1000.;
  logminmass = log(minmass);
  logmaxmass = log(maxmass);
  n_bin = 50;

  c_massfunction = calloc(n_bin,sizeof(MyIDtype));
  massfunction = calloc(n_bin,sizeof(double));
  m= calloc(n_bin,sizeof(double));
  logbin =  (logmaxmass - logminmass)/ n_bin;
 
  for(ihalo = 0; ihalo < TotNhalos; ihalo++)
    {
      if(SubTree[ihalo] == NULLPOINT)
	{
	  Mass = HaloTable[ihalo].Vmax;
	  if(Mass > 1.e-30)
	    {
	      block = (int) ((log(Mass)-logminmass)/logbin);
	      if( block < n_bin)
		{
		  c_massfunction[block]++;
		}
	    }
	}
    }
  usedhalos = 0;
  for(i=0;i<n_bin;i++)
    {
      usedhalos += c_massfunction[i];
    }
  for(i=0;i<n_bin;i++)
    {
      m[i] = exp(logbin*i + logminmass);
      massfunction[i] = (double) c_massfunction[i]/ (double) usedhalos;
      //printf("%lf\t%lf\t%llu\n",m[i],massfunction[i],c_massfunction[i]);
    }
  free(m);
  free(c_massfunction);
  free(massfunction);
}
/*---------------------------------------------------------------------
  Analyse properties of mergers
  -----------------------------------------------------------------------
*/
void merger_analysis(float minmass, float maxmass, int highlim_npart)
{
  unsigned int isnap, jsnap, nStep, block;
  MyIDtype i,j,ihalo,jhalo,MaxNmergers;
  double  latemass;
  double binsize;
  long double summass,dummydouble;
  FILE *fp1, *fp2, *fp3;
  char lostmassOutput_host[MAXSTRING],lostmassOutput[MAXSTRING], nMergerOutput[MAXSTRING],dummyfile[MAXSTRING];

  sprintf(lostmassOutput,"%s_%dLostmass.dat",inputFile,highlim_npart);
  sprintf(lostmassOutput_host,"%s_%dLostmass_host.dat",inputFile,highlim_npart);
  sprintf(nMergerOutput,"%s_%dNmergers.dat",inputFile,highlim_npart);

  /* set up bins for histogram */
  binsize = 0.02;
  nStep = (int) (2./binsize);

  MaxNmergers = 30;
  /* initialise output data  */
  printf("allocate nmerger\n");
  nMergers = calloc(MaxNmergers, sizeof(MyIDtype));

  printf("allocate pp\n");

  pp_nmerger = calloc(MaxNmergers,sizeof(MyIDtype*));
  for(i=0;i<MaxNmergers;i++)
    {
      pp_nmerger[i] = calloc(1,sizeof(MyIDtype));
      pp_nmerger[i][0] = NULLPOINT;
    }
  lostmass = calloc(nStep, sizeof(MyIDtype));
#ifdef READPARTICLE
  if((output.outputFormat > 0.999 && output.outputFormat < 1.001) || (output.outputFormat > 1.129 && output.outputFormat < 1.131))
    {

      lostmass_host = calloc(nStep, sizeof(MyIDtype));
      pp_lostmass = calloc(nStep,sizeof(MyIDtype*));
      for(i=0;i<nStep;i++)
	{
	  pp_lostmass[i] = calloc(1,sizeof(MyIDtype));
	  pp_lostmass[i][0] = NULLPOINT;
	}
    }
#endif
  /* Loop over all halos in the last snapshot  */
  printf("start recursive\n");
  for(i=TotNhalos-SnapNhalos[LASTSNAP]; i < TotNhalos; i++)
    {
      ihalo = i;
#ifdef VERBOSE2
      printf("ihalo = %llu\n",ihalo);
#endif
      if(HaloTable[ihalo].Mvir > 1.e-30 && HaloTable[ihalo].npart > 0)
	{
	  count_mergers(ihalo,nStep,binsize,minmass,maxmass);
	}
#ifdef VERBOSE2
      printf("finish ihalo = %llu\n",ihalo);
#endif
    }
  printf("stop recursive\n");
  fp1 = fopen(lostmassOutput, "w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%f\t%llu\n",((float)(i)+.5)*binsize/2.,lostmass[i]);
    }
  fclose(fp1);
#ifdef READPARTICLE
  if((output.outputFormat > 0.999 && output.outputFormat < 1.001) || (output.outputFormat > 1.129 && output.outputFormat < 1.131))
    {
      fp2 = fopen(lostmassOutput_host, "w+");
      for(i=0;i<nStep;i++)
	{
	  fprintf(fp2,"%f\t%llu\n",((float)(i)+.5)*binsize,lostmass_host[i]);
	}
      fclose(fp2);

      sprintf(dummyfile,"%s_%dLostmass_bin.dat",inputFile,highlim_npart);
      fp1 = fopen(dummyfile,"w+");
      for(i=0;i<nStep;i++)
	{
	  fprintf(fp1,"%f\t",((float)(i)+.5)*binsize);
	  for(j=0;j<lostmass_host[i];j++)
	    {
	      if(pp_lostmass[i][j] < NULLPOINT)
		fprintf(fp1,"%llu\t",pp_lostmass[i][j]);
	    }
	  fprintf(fp1,"\n");
	  free(pp_lostmass[i]);
	}
      fclose(fp1);
      free(pp_lostmass);
      free(lostmass_host);
    }
#endif
  printf("print out nmerger");

  sprintf(dummyfile,"%s_%dNmergers_bin.dat",inputFile,highlim_npart);
  fp1 = fopen(dummyfile,"w+");
  for(i=0;i<MaxNmergers;i++)
    {
      fprintf(fp1,"%d\t",i);
      for(j=0;j<nMergers[i];j++)
	{
	  if(pp_nmerger[i][j] < NULLPOINT)
	    fprintf(fp1,"%llu\t",pp_nmerger[i][j]);
	}
      fprintf(fp1,"\n");
      free(pp_nmerger[i]);
    }
  fclose(fp1);
  free(pp_nmerger);
  fp1 = fopen(nMergerOutput, "w+");
  for(i=0;i<MaxNmergers;i++)
    {
      fprintf(fp1,"%llu\t%llu\n",i,nMergers[i]);
    }
  fclose(fp1);
  free(nMergers);
  free(lostmass);
}

void count_mergers(MyIDtype haloid, int nStep, double binsize,float minmass, float maxmass)
{
  MyIDtype j,k,ihalo_loop,count,jhalo_loop,ihalo,jhalo,khalo,count_particle, temp_count_particle;
  MyIDtype *inputlist,*outputlist,ref_uid,gainpart,lostpart,jid,kid;
  double mass_in, mass_out, mass_lost;
  int block,i;
  unsigned int isnap, jsnap;
  
  //printf("counting mergers for %llu\n",haloid);

  if(output.progs[haloid].nProgs > 0 )
    {
      count = 0;
      mass_out = HaloTable[haloid].Mvir;
      isnap = HaloTable[haloid].SnapID;
      jhalo = output.progs[haloid].progID[0];
      jsnap = HaloTable[jhalo].SnapID;
      mass_in = 0.;
#ifdef READPARTICLE  
      if((output.outputFormat > 0.999 && output.outputFormat < 1.001) || (output.outputFormat > 1.129 && output.outputFormat < 1.131))
	{
#ifdef VERBOSE2
	  printf("malloc inputlist\n");
#endif
	  inputlist = malloc(0);
#ifdef VERBOSE2
	  printf("malloc inputlist completed\n");
#endif
	  count_particle = 0;
	}
#endif
      for(ihalo_loop = 0; ihalo_loop < output.progs[haloid].nProgs; ihalo_loop++)
	{
	  jhalo_loop = output.progs[haloid].progID[ihalo_loop];
#ifdef VERBOSE2
	  printf("jhalo_loopX %llu : %llu particles\n",jhalo_loop,HaloTable[jhalo_loop].npart);
#endif
	  if(HaloTable[jhalo_loop].Mvir > 1.e-30 && HaloTable[jhalo_loop].npart > 0)
	    {

	      count_mergers(jhalo_loop,nStep,binsize,minmass,maxmass);
	      mass_in += HaloTable[jhalo_loop].Mvir;
	      count++;
#ifdef READPARTICLE
	      if((output.outputFormat > 0.999 && output.outputFormat < 1.001) || (output.outputFormat > 1.129 && output.outputFormat < 1.131))
		{
#ifdef VERBOSE2
		  printf("halo %llu : %llu particles\n",jhalo_loop,HaloTable[jhalo_loop].npart);
#endif
		  count_particle += HaloTable[jhalo_loop].npart;
#ifdef VERBOSE2
		  printf("realloc inputlist\n");
#endif
		  inputlist = realloc(inputlist,count_particle*sizeof(MyIDtype));
#ifdef VERBOSE2
		  printf("realloc inputlist completed\n");
#endif
		  for(j=count_particle-HaloTable[jhalo_loop].npart;j<count_particle;j++)
		    {
		      inputlist[j] = HaloTable[jhalo_loop].Particles[j].ParticleID;
		    }
		}
#endif
	    }
	}

      nMergers[count]++;
      //printf("realloc p_nmerger %d\n",count);
      pp_nmerger[count] = realloc(pp_nmerger[count],nMergers[count]*sizeof(MyIDtype));
      pp_nmerger[count][nMergers[count]-1] = haloid;
      // if they use original catalogue
#ifdef READPARTICLE
#ifdef VERBOSE2
      printf("start looping\n");
#endif
      if((output.outputFormat > 0.999 && output.outputFormat < 1.001) || (output.outputFormat > 1.129 && output.outputFormat < 1.131))
	{
	  qsort(inputlist,count_particle, sizeof(MyIDtype), compareMyIDType);
	  ref_uid = inputlist[0];
	  for(j=1;j<count_particle;j++)
	    {
	      if(inputlist[j] == ref_uid)
		{
		  inputlist[j] = NULLPOINT;
		}
	      else
		{
		  ref_uid = inputlist[j];
		}
	    }
	  qsort(inputlist,count_particle, sizeof(MyIDtype), compareMyIDType);
	  temp_count_particle = 0;
	  for(j=0;j<count_particle;j++)
	    {
	      if(inputlist[j] < NULLPOINT)
		{
		  temp_count_particle++;
		}
	      else
		break;
	    }
	  count_particle = temp_count_particle;
#ifdef VERBOSE2
	  printf("realloc in/outputlist\n");
#endif
	  inputlist = realloc(inputlist,temp_count_particle*sizeof(MyIDtype));
	  outputlist = malloc(HaloTable[haloid].npart*sizeof(MyIDtype));
#ifdef VERBOSE2
	  printf("finish realloc in/outputlist\n");
#endif
	  for(j=0;j<HaloTable[haloid].npart;j++)
	    {
	      outputlist[j] = HaloTable[haloid].Particles[j].ParticleID;
	    }
	  qsort(outputlist,HaloTable[haloid].npart, sizeof(MyIDtype), compareMyIDType);
	  // very long winded to compare particle list
	  gainpart = 0;
	  lostpart = 0;
#ifdef VERBOSE2
	  printf("start omp loop\n");
#endif
	  if(SubTree[haloid] == NULLPOINT && SubTree[jhalo] == NULLPOINT)
	    {
	      if(HaloTable[jhalo].Mvir > minmass && HaloTable[jhalo].Mvir < maxmass && HaloTable[haloid].Mvir > minmass && HaloTable[haloid].Mvir < maxmass)
		{
#pragma omp parallel for default(shared) private(jid) reduction(+:lostpart) 
		  for(j = 0; j< count_particle; j++)
		    {
		      jid = inputlist[j];
		      if(Generalsearch(jid,HaloTable[haloid].npart,outputlist) < NULLPOINT)
			{
			  lostpart++;
			  //printf("haloid = %llu\n",haloid);
			}
		      /*
			for(k = 0; k < HaloTable[haloid].npart; k++)
			{
			kid = outputlist[k];
			if(kid == jid)
			{
			lostpart++;
			break;
			}
			}
		      */
		    }
		}
	    }
#ifdef VERBOSE2
	  printf("stop omp loop\n");
#endif
	  lostpart = count_particle - lostpart;
	  /*
	    for(k = 0; k <HaloTable[haloid].npart; k++)
	    {
	    kid = outputlist[k];
	    //printf("kid = %llu\n",kid);
	    for(j = 0; j< count_particle; j++)
	    {
	    jid = inputlist[j];
	    if(kid == jid)
	    {
	    gainpart++;
	    break;
	    }
	    }
	    }
	    gainpart = HaloTable[haloid].npart-gainpart;
	  */
	  free(inputlist);
	  free(outputlist);
#ifdef VERBOSE2
	  printf("finish free\n");
#endif
	  //printf("lostpart = %llu/%llu \n",lostpart,count_particle);
	  //exit(0);
	  //mass_lost =  (double)lostpart /(double)(count_particle + HaloTable[haloid].npart)/(snapTime[isnap] - snapTime[jsnap])*(snapTime[isnap] + snapTime[jsnap]);
#ifdef VERBOSE2
	  printf("lostpart = %llu, count_particle = %llu\n", lostpart,count_particle);
#endif
	  mass_lost =  (double)lostpart /(double)(count_particle )*2.;


	  if(count_particle > 0)
	    {
	      //block = (int) ((atan(mass_lost)/asin(1.) +1.)/binsize);
	      block = (int) ( mass_lost  /  binsize);
	      if(block == nStep)
		block = nStep-1;
#ifdef VERBOSE2	  
	      printf("record block = %d\n",block);
#endif
	      if(SubTree[haloid] == NULLPOINT && SubTree[jhalo] == NULLPOINT)
		{
		  if(HaloTable[jhalo].Mvir > minmass && HaloTable[jhalo].Mvir < maxmass && HaloTable[haloid].Mvir > minmass && HaloTable[haloid].Mvir < maxmass)
		    {
		      lostmass_host[block]++;
		      //printf("realloc pp_lostmass %d\n",block);
		      pp_lostmass[block] = realloc(pp_lostmass[block],lostmass_host[block]*sizeof(MyIDtype));
		      pp_lostmass[block][lostmass_host[block]-1] = haloid;
		    }
		}  
	    }
	}
#ifdef VERBOSE2
      printf("complete one loop\n");
#endif

#endif // READPARTICLE
      mass_lost = (mass_out - mass_in)/(mass_out+mass_in)/(snapTime[isnap] - snapTime[jsnap])*(snapTime[isnap] + snapTime[jsnap]);
      block = (int) ((atan(mass_lost)/asin(1.) +1.)/binsize);
      if(SubTree[haloid] == NULLPOINT && SubTree[jhalo] == NULLPOINT)
	{
	  if(HaloTable[jhalo].Mvir > minmass && HaloTable[jhalo].Mvir < maxmass && HaloTable[haloid].Mvir > minmass && HaloTable[haloid].Mvir < maxmass)
	    lostmass[block]++;
	} 
    }
  else
    {
      count = 0;
      if(HaloTable[haloid].SnapID > 30)
	{
	  nMergers[count]++;
	  //printf("realloc p_nmerger %d\n",count);
	  pp_nmerger[count] = realloc(pp_nmerger[count],nMergers[count]*sizeof(MyIDtype));
	  //printf("set haloid=%llu %d %llu\n",haloid,count,nMergers[count]-1);
	  pp_nmerger[count][nMergers[count]-1] = haloid;
	}
      //printf("finish loop\n");
    }
}
/*---------------------------------------------------------------------
  Analyse properties of main branches between snapshots
  -----------------------------------------------------------------------
*/
void snapshot_stats(int snapid, float minmass, float maxmass, int highlim_npart)
{
  unsigned int count_id,nStep,iSnap,jSnap,kSnap,block,count,savedblock;
  MyIDtype i,j,k,ihalo,jhalo,khalo,sumhalo,start_hid;
  double binsize,lastmass,corr_self,firstslope,secondslope,saveddouble,massboundery,Vmax,distance;
  double V_j[3],V_i[3],R_i[3],R_j[3],Rvir,Rmax,disp[3],expected_R[3],avg_V[3],sum_R;
  long double Radius;
  MyIDtype *MeanBeta;
  MyIDtype *MeanBeta_l,*MeanBeta_h,*MeanBeta_host;
  MyIDtype **pp_MeanBeta;
  MyIDtype **MeanBeta_l_pp,**MeanBeta_h_pp,**MeanBeta_host_pp;
  MyIDtype *dispCorr_host;
  MyIDtype **dispCorr_host_pp;
  MyIDtype *CorrBeta_h;
  MyIDtype **CorrBeta_h_pp;
  MyIDtype *Corr_h;
  MyIDtype **Corr_h_pp;
  MyIDtype *nmergers;
  char meanbetaOutput[MAXSTRING];
  char l_meanbetaOutput[MAXSTRING],h_meanbetaOutput[MAXSTRING],host_meanbetaOutput[MAXSTRING],h_meanbetaOutput_pp[MAXSTRING];
  char host_dispOutput[MAXSTRING],host_dispOutput_pp[MAXSTRING];
  char betacorr_highmass[MAXSTRING], betacorr_highmass_pp[MAXSTRING];
  char corr_highmass[MAXSTRING], corr_highmass_pp[MAXSTRING];
  char nmergers_output[MAXSTRING];
  int *flag, checkref;
  FILE *rc1;
  FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8, *fp9;
  printf("Start analysis for snapshot %d\n",snapid);
  /*
    Set up filenames for outputs
  */

  // set highmass/ low mass at 100 particle masses
  massboundery = p_mass*(float) highlim_npart;

  // set binsize for histograms
  binsize = 0.02;
  nStep = (int) (2./binsize);

  start_hid = 0;
  for(count_id=0;count_id<snapid;count_id++)
    {
      start_hid += SnapNhalos[count_id];
    }

  MeanBeta = calloc(nStep,sizeof(MyIDtype));
  MeanBeta_l = calloc(nStep,sizeof(MyIDtype));
  MeanBeta_h = calloc(nStep,sizeof(MyIDtype));
  MeanBeta_host = calloc(nStep,sizeof(MyIDtype));
  dispCorr_host = calloc(nStep,sizeof(MyIDtype));
  CorrBeta_h = calloc(nStep,sizeof(MyIDtype));
  Corr_h = calloc(nStep,sizeof(MyIDtype));
  nmergers= calloc(30,sizeof(MyIDtype));

  MeanBeta_h_pp = calloc(nStep,sizeof(MyIDtype *));
  dispCorr_host_pp = calloc(nStep,sizeof(MyIDtype *));
  CorrBeta_h_pp = calloc(nStep,sizeof(MyIDtype *));
  Corr_h_pp = calloc(nStep,sizeof(MyIDtype *));

  for(i=0;i<nStep;i++)
    {
      dispCorr_host_pp[i] = calloc(1,sizeof(MyIDtype));
      dispCorr_host_pp[i][0] = NULLPOINT;
      CorrBeta_h_pp[i] = calloc(1,sizeof(MyIDtype));
      CorrBeta_h_pp[i][0] = NULLPOINT;
      Corr_h_pp[i] = calloc(1,sizeof(MyIDtype));
      Corr_h_pp[i][0] = NULLPOINT;
      MeanBeta_h_pp[i] = calloc(1,sizeof(MyIDtype));
      MeanBeta_h_pp[i][0] = NULLPOINT;
    }
  printf("Finish allocate\n");
  /*
    Loop over halos in the last snapshot
  */
  for(i=start_hid; i < start_hid+SnapNhalos[snapid]; i++)
    {
      //printf("i = %llu\n",i);
      /* record properties of the halo in the last snapshot */
      lastmass = HaloTable[i].Mvir; 
      if(lastmass > 1.e-35)
	{
	  ihalo = i;
	  iSnap = HaloTable[ihalo].SnapID; 
	  count = 0;
	  for(k=0;k<output.progs[ihalo].nProgs;k++)
	    {
	      jhalo = output.progs[ihalo].progID[k];
	      if(HaloTable[jhalo].Mvir > 1.e-35)
		{
		  count++;
		}
	    }
	  nmergers[count]++;

	  if(output.progs[ihalo].nProgs > 0 )
	    {
	      /* snapshot ID for ihalo  */
	      iSnap = HaloTable[ihalo].SnapID;

	      /* set the main progenitor of ihalo to be jhalo */
	      jhalo = output.progs[ihalo].progID[0];
	  
	      jSnap = HaloTable[jhalo].SnapID;
	  
	      if(HaloTable[jhalo].Mvir > 1.e-35)
		{
		  /* 6D space for ihalo */
		  R_i[0] = HaloTable[ihalo].Xc;
		  R_i[1] = HaloTable[ihalo].Yc;
		  R_i[2] = HaloTable[ihalo].Zc;

		  V_i[0] = HaloTable[ihalo].VXc;
		  V_i[1] = HaloTable[ihalo].VYc;
		  V_i[2] = HaloTable[ihalo].VZc;


		  /* 6D space for jhalo  */
		  R_j[0] = HaloTable[jhalo].Xc;
		  R_j[1] = HaloTable[jhalo].Yc;
		  R_j[2] = HaloTable[jhalo].Zc;

		  V_j[0] = HaloTable[jhalo].VXc;
		  V_j[1] = HaloTable[jhalo].VYc;
		  V_j[2] = HaloTable[jhalo].VZc;

		  /* record Rmin and Rmax for jhalo  */
		  
		  /*
		    calculate delta_r = (r_n - (r_{n-1} +vt)) / Rvir *in units km,s*
		  */

		  Radius = 0.;
		  for(k=0;k<3;k++)
		    {
		      // change a\dot{x} => \dot{x}/h as well (r is in kpc/h)
		      //avg_V[k] = 0.5/h*(V_i[k]/expansion_factor[iSnap]+V_j[k]/expansion_factor[jSnap]);
		      avg_V[k] = V_j[k]/expansion_factor[jSnap];

		      //avg_V[k] = snap_stats[jsnap].Vrms/expansion_factor[jSnap]/h;
		      expected_R[k] = fmod((R_j[k] + avg_V[k]*(snapTimeYear[iSnap] - snapTimeYear[jSnap])*year2second/kpc2km) + Boxsize,Boxsize);
		      sum_R = (avg_V[k])*(snapTimeYear[iSnap] - snapTimeYear[jSnap])*year2second/kpc2km;
		      
		      Radius += (long double) (sum_R * sum_R);
		      //disp[k] = MIN(fabs(R_j[k] - R_i[k]),Boxsize - fabs(R_j[k] - R_i[k]));
		      //disp[k] = MIN(fabs(R_j[k] - R_i[k]),Boxsize - fabs(R_j[k] - R_i[k])) ;
		      disp[k] = MIN( fabs(R_i[k]*kpc2km - expected_R[k]),Boxsize*kpc2km - fabs(R_i[k]*kpc2km - expected_R[k]) ); 
		      //disp[k] /= (Rvir*kpc2km);
		    }
		  //Vmax = MAX(snap_stats[iSnap].Vmax/expansion_factor[iSnap],snap_stats[jSnap].Vmax/expansion_factor[iSnap])/h;
		  Rmax =  sqrt(HaloTable[ihalo].Rmax*HaloTable[ihalo].Rmax + HaloTable[jhalo].Rmax*HaloTable[jhalo].Rmax);
		  //Rvir = sqrt(HaloTable[ihalo].Rvir* HaloTable[ihalo].Rvir + HaloTable[jhalo].Rvir*HaloTable[jhalo].Rvir);
		  Rvir = (HaloTable[ihalo].Rvir+HaloTable[jhalo].Rvir)/2.;
		  Radius = (sqrt(Radius)+Rvir);
		  distance =  sqrt(disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2]);
		  distance /= (double) Radius;
		  //distance /= (Vmax*(snapTimeYear[iSnap] - snapTimeYear[jSnap])*year2second);
	      	  /*
		    calculate arctan(delta_r)/pi
		  */
		  firstslope = atan(distance)/(asin(1.));
		  /*
		    put in histogram
		  */
		  block = (int) ((firstslope+1.)/binsize);
		  //printf("Block %d\n",block);
		  if(block == nStep) 
		    block = nStep-1;

		  if(SubTree[ihalo] == NULLPOINT && SubTree[jhalo] == NULLPOINT)
		    {
		      dispCorr_host[block]++;
		      dispCorr_host_pp[block] = realloc(dispCorr_host_pp[block],dispCorr_host[block]*sizeof(MyIDtype));
		      dispCorr_host_pp[block][dispCorr_host[block]-1] = ihalo;
		    }

		  /* calculate beta */
		  firstslope = atan((HaloTable[ihalo].Mvir - HaloTable[jhalo].Mvir) /(HaloTable[ihalo].Mvir + HaloTable[jhalo].Mvir)/(snapTime[iSnap] - snapTime[jSnap]) * (snapTime[iSnap] + snapTime[jSnap]) )/(asin(1.)); 
		  /* put it histogram */
		  block = (int) ((firstslope+1.)/binsize);
		  //printf("Block %d\n",block);
		  if(block == nStep) 
		    block = nStep-1;
		  /* discard when beta = 0 */
		  MeanBeta[block]++;
		  //printf("block %d N= %d\n",block, MeanBeta[block]);

		  

		  //printf("finish adding element\n");
		  if(HaloTable[ihalo].Mvir < massboundery)
		    MeanBeta_l[block]++;
		  else
		    {
		      MeanBeta_h[block]++;
		      MeanBeta_h_pp[block] = realloc(MeanBeta_h_pp[block],MeanBeta_h[block]*sizeof(MyIDtype));
		      MeanBeta_h_pp[block][MeanBeta_h[block]-1] = ihalo;
		    }
		  //printf("plost L/H\n");
		  if(SubTree[ihalo] == NULLPOINT && SubTree[jhalo] == NULLPOINT)
		    MeanBeta_host[block]++;
		  //printf("plost host\n");
		  //printf("start k\n");
		  firstslope = 2.*(HaloTable[ihalo].Mvir - HaloTable[jhalo].Mvir)/(HaloTable[ihalo].Mvir + HaloTable[jhalo].Mvir);
		  
		  if(output.progs[jhalo].nProgs > 0 )
		    {
		      khalo = output.progs[jhalo].progID[0];
		      kSnap = HaloTable[khalo].SnapID; 
		      if(HaloTable[khalo].Mvir > 1.e-30)
			{
			  firstslope = atan((HaloTable[ihalo].Mvir - HaloTable[khalo].Mvir) / (HaloTable[ihalo].Mvir + HaloTable[khalo].Mvir) /(snapTime[iSnap] - snapTime[kSnap]) * (snapTime[iSnap] + snapTime[kSnap]) ); 
			  /* calculate beta_{n-1,n} */
			  secondslope = atan((HaloTable[jhalo].Mvir - HaloTable[khalo].Mvir)/ (HaloTable[jhalo].Mvir + HaloTable[khalo].Mvir) /(snapTime[jSnap] - snapTime[kSnap]) * (snapTime[jSnap] + snapTime[kSnap]));

			  /* beta_{n-1,n} - beta_{n-1,n+1} */
			  corr_self = fabs(secondslope-firstslope)/asin(1.)/2.;
	      
			  //printf("corr %lf\n",corr_self);
			  block = (int) ((corr_self + 1.)/binsize);
			  //printf("Block %d\n",block);
			  if(block == nStep) 
			    block = nStep-1;

			  if(HaloTable[ihalo].Mvir > massboundery && HaloTable[jhalo].Mvir > massboundery && HaloTable[khalo].Mvir > massboundery)
			    {			
			      CorrBeta_h[block]++;
			      CorrBeta_h_pp[block] = realloc(CorrBeta_h_pp[block],CorrBeta_h[block]*sizeof(MyIDtype));
			      CorrBeta_h_pp[block][CorrBeta_h[block]-1] = ihalo;
			    }
			  firstslope = atan((HaloTable[ihalo].Mvir - HaloTable[khalo].Mvir) / (HaloTable[ihalo].Mvir + HaloTable[khalo].Mvir) /(snapTime[iSnap] - snapTime[kSnap]) * 2.); 
			  /* calculate beta_{n-1,n} */
			  secondslope = atan((HaloTable[jhalo].Mvir - HaloTable[khalo].Mvir)/ (HaloTable[jhalo].Mvir + HaloTable[khalo].Mvir) /(snapTime[jSnap] - snapTime[kSnap]) * 2.);

			  /* beta_{n-1,n} - beta_{n-1,n+1} */
			  corr_self = fabs(secondslope-firstslope)/asin(1.)/2.;
	      
			  //printf("corr %lf\n",corr_self);
			  block = (int) ((corr_self + 1.)/binsize);
			  //printf("Block %d\n",block);
			  if(block == nStep) 
			    block = nStep-1;
			   if(HaloTable[ihalo].Mvir > massboundery && HaloTable[jhalo].Mvir > massboundery && HaloTable[khalo].Mvir > massboundery)
			    {			
			      Corr_h[block]++;
			      Corr_h_pp[block] = realloc(Corr_h_pp[block],Corr_h[block]*sizeof(MyIDtype));
			      Corr_h_pp[block][Corr_h[block]-1] = ihalo;
			    }

			}
		    }
		  //printf("start k\n");
		}
	    }
	}
    }

  printf("Set up filenames\n");
  sprintf(meanbetaOutput,"%s_%dMeanBeta.%d.dat",inputFile,highlim_npart,snapid);
  sprintf(l_meanbetaOutput,"%s_%dMeanBetaLowmass.%d.dat",inputFile,highlim_npart,snapid);
  sprintf(h_meanbetaOutput,"%s_%dMeanBetaHighmass.%d.dat",inputFile,highlim_npart,snapid);
  sprintf(host_meanbetaOutput,"%s_%dMeanBetaHost.%d.dat",inputFile,highlim_npart,snapid);
  sprintf(host_dispOutput,"%s_%dDispMainHost.%d.dat",inputFile,highlim_npart,snapid);
  sprintf(nmergers_output,"%s_%dNmergers.%d.dat",inputFile,highlim_npart,snapid);

  sprintf(betacorr_highmass,"%s_%dBetaCorrHighmass.%d.dat",inputFile,highlim_npart,snapid );
  sprintf(corr_highmass,"%s_%dSlopeCorrHighmass.%d.dat",inputFile,highlim_npart,snapid);

  sprintf(host_dispOutput_pp,"%s/%s_%dDispMainHost_pp.%d.dat",ppFolder,inputFile,highlim_npart,snapid);
  sprintf(betacorr_highmass_pp,"%s/%s_%dBetaCorrHighmass_pp.%d.dat",ppFolder,inputFile,highlim_npart,snapid);
  sprintf(h_meanbetaOutput_pp,"%s/%s_%dMeanBetaHighmass_pp.%d.dat",ppFolder,inputFile,highlim_npart,snapid);
  sprintf(corr_highmass_pp,"%s/%s_%dSlopeCorrHighmass.%d.dat",ppFolder,inputFile,highlim_npart,snapid);

  printf("finish make new file names\n");

  fp1 = fopen(nmergers_output,"w+");

  for(i=0;i<30;i++)
    {
      fprintf(fp1,"%d\t%llu\n",i,nmergers[i]);
    }
  fclose(fp1);

  fp2 = fopen(meanbetaOutput,"w+");
  fp3 = fopen(host_dispOutput,"w+");

  for(i=0;i<nStep;i++)
    {
      fprintf(fp2,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,MeanBeta[i]);
      fprintf(fp3,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,dispCorr_host[i]);
    }
  fclose(fp2);
  fclose(fp3);



  fp1 = fopen(host_dispOutput_pp,"w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%d\t", (int) i);
      for(j=0;j<dispCorr_host[i];j++)
	{
	  if( dispCorr_host_pp[i][j] < NULLPOINT ) 
	    fprintf(fp1,"%llu\t",dispCorr_host_pp[i][j]);
	}
      fprintf(fp1,"\n");
      free(dispCorr_host_pp[i]);
    }
  fclose(fp1);

  fp1 = fopen(betacorr_highmass_pp,"w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%d\t", (int) i);
      for(j=0;j <CorrBeta_h[i]; j++)
	{
	  if(CorrBeta_h_pp[i][j] < NULLPOINT ) 
	    fprintf(fp1,"%llu\t",CorrBeta_h_pp[i][j]);
	}
      fprintf(fp1,"\n");
      free(CorrBeta_h_pp[i]);
    }
  fclose(fp1);

  fp1 = fopen(corr_highmass_pp,"w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%d\t", (int) i);
      for(j=0;j <Corr_h[i]; j++)
	{
	  if(Corr_h_pp[i][j] < NULLPOINT ) 
	    fprintf(fp1,"%llu\t",Corr_h_pp[i][j]);
	}
      fprintf(fp1,"\n");
      free(Corr_h_pp[i]);
    }
  fclose(fp1);
    
  fp1 = fopen(h_meanbetaOutput_pp,"w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%d\t", (int) i);
      for(j=0;j < MeanBeta_h[i]; j++)
	{
	  if(MeanBeta_h_pp[i][j] < NULLPOINT ) 
	    fprintf(fp1,"%llu\t",MeanBeta_h_pp[i][j]);
	}
      fprintf(fp1,"\n");
      free(MeanBeta_h_pp[i]);
    }
  fclose(fp1);

  printf("finish lot 1\n");
  fp1 = fopen(l_meanbetaOutput,"w+");
  fp2 = fopen(h_meanbetaOutput,"w+");
  fp3 = fopen(host_meanbetaOutput, "w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,MeanBeta_l[i]);
      fprintf(fp2,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,MeanBeta_h[i]);
      fprintf(fp3,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,MeanBeta_host[i]);
    }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);

  fp1 = fopen(betacorr_highmass,"w+");
  fp2 = fopen(corr_highmass,"w+");
 
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,CorrBeta_h[i]);
      fprintf(fp2,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,Corr_h[i]);

    }
  fclose(fp1);
  fclose(fp2);


  free(MeanBeta);
  free(MeanBeta_l);
  free(MeanBeta_host);

  free(dispCorr_host);
  free(dispCorr_host_pp);
  free(CorrBeta_h);
  free(CorrBeta_h_pp);
  free(Corr_h);
  free(Corr_h_pp);
  free(MeanBeta_h);
  free(MeanBeta_h_pp);
  printf("finish lot 2\n");
}

void main_branch_analysis(float minmass, float maxmass, int highlim_npart)
{
  unsigned int snapid,nStep,iSnap,jSnap,kSnap,block,count,savedblock;
  MyIDtype i,j,k,ihalo,jhalo,khalo,sumhalo;
  double binsize,lastmass,corr_self,firstslope,secondslope,saveddouble,massboundery,Vmax,distance;
  double V_j[3],V_i[3],R_i[3],R_j[3],R_k[3],Rvir,Rmax,disp[3],expected_R[3],avg_V[3],sum_R,R1[3],R2[3],abs_R1,abs_R2;
  long double Radius;
  MyIDtype *Corr,*CorrBeta,*Mean,*MeanBeta,*linkdepth,*dispCorr, *AlphaSub, *BetaSub,*highmasslinkdepth,*lowmasslinkdepth,*hostlinkdepth,*resilinkdepth;
  MyIDtype *Mean_l,*Mean_h;
  MyIDtype *MeanBeta_l,*MeanBeta_h;
  MyIDtype *Corr_l,*Corr_h;
  MyIDtype *CorrBeta_l,*CorrBeta_h;
  MyIDtype *dispCorr_l, *dispCorr_h;
  MyIDtype *dispCorr_resi, *dispCorr_host;
  MyIDtype *postFluct;
  MyIDtype **linkdepth_pp;
  MyIDtype **beta_pp;
  char meanOutput[MAXSTRING],corrOutput[MAXSTRING],meanbetaOutput_pp[MAXSTRING],meanbetaOutput[MAXSTRING],depthOutput[MAXSTRING],mainbranchmassOutput[MAXSTRING],dispOutput[MAXSTRING];
  char corrbetaOutput[MAXSTRING],corrbetaOutput_pp[MAXSTRING],alphasubOutput[MAXSTRING],betasubOutput[MAXSTRING];
  char lowmassLinkdepth[MAXSTRING],highmassLinkdepth[MAXSTRING],host_Linkdepth[MAXSTRING],resi_Linkdepth[MAXSTRING];
  char l_meanOutput[MAXSTRING],h_meanOutput[MAXSTRING];
  char l_meanbetaOutput[MAXSTRING],h_meanbetaOutput[MAXSTRING];
  char l_corrOutput[MAXSTRING],h_corrOutput[MAXSTRING];
  char l_corrbetaOutput[MAXSTRING],h_corrbetaOutput[MAXSTRING];
  char l_depthOutput[MAXSTRING], h_depthOutput[MAXSTRING];
  char l_dispOutput[MAXSTRING],h_dispOutput[MAXSTRING];
  char resi_dispOutput[MAXSTRING],host_dispOutput[MAXSTRING],host_dispOutput_pp[MAXSTRING];
  char postFluctOutput[MAXSTRING];
  char linkdepthOutput_pp[MAXSTRING];
  double* meanmainbranch;
  MyIDtype* TotMainbranch;
  
  int *flag, checkref;
  FILE *rc1;
  FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8, *fp9;
  printf("Start main tree analysis\n");
  /*
    Set up filenames for outputs
  */
  printf("Start setting filenames\n");

  sprintf(host_dispOutput_pp,"%s_%dDispMainHost_pp.dat",inputFile,highlim_npart);
  sprintf(corrbetaOutput_pp,"%s_%dBeta_pp.dat",inputFile,highlim_npart);
  sprintf(linkdepthOutput_pp,"%s_%dLinkdepth_pp.dat",inputFile,highlim_npart);

  // set highmass/ low mass at 100 particle masses
  massboundery = p_mass*(float) highlim_npart;

  // set binsize for histograms
  binsize = 0.02;
  nStep = (int) (2./binsize);

  /*
    set up initial values for output data
  */

  Corr = calloc(nStep,sizeof(MyIDtype));
  Corr_l = calloc(nStep,sizeof(MyIDtype));
  Corr_h = calloc(nStep,sizeof(MyIDtype));

  CorrBeta = calloc(nStep,sizeof(MyIDtype));
  CorrBeta_l = calloc(nStep,sizeof(MyIDtype));
  CorrBeta_h = calloc(nStep,sizeof(MyIDtype));

  MeanBeta = calloc(nStep,sizeof(MyIDtype));
  MeanBeta_l = calloc(nStep,sizeof(MyIDtype));
  MeanBeta_h = calloc(nStep,sizeof(MyIDtype));

  Mean = calloc(nStep,sizeof(MyIDtype));
  Mean_l = calloc(nStep,sizeof(MyIDtype));
  Mean_h = calloc(nStep,sizeof(MyIDtype));

  linkdepth = calloc(NSNAPS,sizeof(MyIDtype));

  dispCorr = calloc(nStep,sizeof(MyIDtype));
  dispCorr_l = calloc(nStep,sizeof(MyIDtype));
  dispCorr_h = calloc(nStep,sizeof(MyIDtype));

  dispCorr_resi = calloc(nStep,sizeof(MyIDtype));
  dispCorr_host = calloc(nStep,sizeof(MyIDtype));
  postFluct = calloc(nStep,sizeof(MyIDtype));

  lowmasslinkdepth = calloc(NSNAPS,sizeof(MyIDtype));
  highmasslinkdepth = calloc(NSNAPS,sizeof(MyIDtype));
  hostlinkdepth = calloc(NSNAPS,sizeof(MyIDtype));
  resilinkdepth = calloc(NSNAPS,sizeof(MyIDtype));

  meanmainbranch = calloc(NSNAPS,sizeof(double));
  TotMainbranch = calloc(NSNAPS,sizeof(MyIDtype));

  //flag = calloc(TotNhalos,sizeof(int));
  //printf("Use nStep:%d\tbinsize:%f\n",nStep,binsize);

  fp2 = fopen(host_dispOutput_pp, "w+");
  fp1 = fopen(linkdepthOutput_pp, "w+");
  fp3 = fopen(corrbetaOutput_pp, "w+");

  /*
    Loop over halos in the last snapshot
  */
  printf("Start loop\n");
  for(i=TotNhalos-SnapNhalos[LASTSNAP]; i < TotNhalos; i++)
    {
      /* record properties of the halo in the last snapshot */
      lastmass = HaloTable[i].Mvir; 
      if(lastmass > 1.e-35)
	{
	  count = 0; // initial the counter for link depth
	  ihalo = i;
	  iSnap = HaloTable[ihalo].SnapID; 
	  if(SubTree[i] == NULLPOINT)
	    {
	      meanmainbranch[iSnap] += HaloTable[ihalo].Mvir/lastmass;
	      TotMainbranch[iSnap]++; // increase number of halo in iSnap 
	    }
	  while(output.progs[ihalo].nProgs > 0 )
	    {
	      /* snapshot ID for ihalo  */
	      iSnap = HaloTable[ihalo].SnapID;

	      /* set the main progenitor of ihalo to be jhalo */
	      jhalo = output.progs[ihalo].progID[0];
	  
	      jSnap = HaloTable[jhalo].SnapID;
	  
	      if(HaloTable[jhalo].Mvir > 1.e-35)
		{
		  /* 6D space for ihalo */
		  R_i[0] = HaloTable[ihalo].Xc;
		  R_i[1] = HaloTable[ihalo].Yc;
		  R_i[2] = HaloTable[ihalo].Zc;

		  V_i[0] = HaloTable[ihalo].VXc;
		  V_i[1] = HaloTable[ihalo].VYc;
		  V_i[2] = HaloTable[ihalo].VZc;

	  
	  
		  count+= iSnap-jSnap; // increase link depth
		  /* increase mass in jsnap  */
		  if(SubTree[jhalo]==NULLPOINT)
		    {
		      meanmainbranch[jSnap] += HaloTable[jhalo].Mvir/lastmass;
		      TotMainbranch[jSnap]++;
		    }

		  /* 6D space for jhalo  */
		  R_j[0] = HaloTable[jhalo].Xc;
		  R_j[1] = HaloTable[jhalo].Yc;
		  R_j[2] = HaloTable[jhalo].Zc;

		  V_j[0] = HaloTable[jhalo].VXc;
		  V_j[1] = HaloTable[jhalo].VYc;
		  V_j[2] = HaloTable[jhalo].VZc;

		  /* record Rmin and Rmax for jhalo  */
		  
	

		  /*
		    calculate delta_r = (r_n - (r_{n-1} +vt)) / Rvir *in units km,s*
		  */

		  Radius = 0.;
		  for(k=0;k<3;k++)
		    {
		      // change a\dot{x} => \dot{x}/h as well (r is in kpc/h)
		      avg_V[k] = 0.5*h*(V_i[k]/expansion_factor[iSnap]+V_j[k]/expansion_factor[jSnap]);
		      //avg_V[k] = V_j[k]/expansion_factor[jSnap];
		      //avg_V[k] = 0.5*h*(V_i[k]+V_j[k]);
		      //avg_V[k] = snap_stats[jsnap].Vrms/expansion_factor[jSnap]/h;
		      expected_R[k] = fmod((R_j[k] + avg_V[k]*(snapTimeYear[iSnap] - snapTimeYear[jSnap])*year2second/kpc2km) + Boxsize,Boxsize);
		      sum_R = (avg_V[k])*(snapTimeYear[iSnap] - snapTimeYear[jSnap])*year2second/kpc2km;
		      //sum_R = fabs(V_i[k]/expansion_factor[iSnap]-V_j[k]/expansion_factor[jSnap])*(snapTimeYear[iSnap] - snapTimeYear[jSnap])*year2second/kpc2km;
		      Radius += (long double) (sum_R * sum_R);
		      //disp[k] = MIN(fabs(R_j[k] - R_i[k]),Boxsize - fabs(R_j[k] - R_i[k]));
		    
		      disp[k] = MIN( fabs(R_i[k] - expected_R[k]),Boxsize - fabs(R_i[k] - expected_R[k]) ); 
		      //disp[k] /= (Rvir*kpc2km);
		    }
		  //Vmax = MAX(snap_stats[iSnap].Vmax/expansion_factor[iSnap],snap_stats[jSnap].Vmax/expansion_factor[iSnap])/h;
		  Rmax =  HaloTable[ihalo].Rmax + HaloTable[jhalo].Rmax;
		  Rvir = (HaloTable[ihalo].Rvir + HaloTable[jhalo].Rvir)/2.;
		  //Radius = 0.;
		  Radius = (sqrt(Radius)+Rvir);
		  distance =  sqrt(disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2]);
		  distance /= (double) Radius;
		  //distance /= (Vmax*(snapTimeYear[iSnap] - snapTimeYear[jSnap])*year2second);
		 

		  /*
		    calculate arctan(delta_r)/pi
		  */
		  firstslope = atan(distance)/(asin(1.));
		  /*
		    put in histogram
		  */
		  block = (int) ((firstslope+1.)/binsize);
		  //printf("Block %d\n",block);
		  if(block == nStep) 
		    block = nStep-1;

		  if(iSnap == 61 ) //this is Yao's case, only the last two snapshots
		    dispCorr[block]++;

		  if(HaloTable[ihalo].Mvir < massboundery)
		    dispCorr_l[block]++;
		  else
		    dispCorr_h[block]++;

		  if(SubTree[ihalo] < NULLPOINT)
		    dispCorr_resi[block]++;
		  
		  if(SubTree[ihalo] == NULLPOINT && SubTree[jhalo] == NULLPOINT)
		    {
		      if(HaloTable[ihalo].Mvir < maxmass && HaloTable[ihalo].Mvir > minmass && HaloTable[jhalo].Mvir < maxmass && HaloTable[jhalo].Mvir > minmass)
			{
			  dispCorr_host[block]++;
			  if(distance > 0.5)
			    {
			      fprintf(fp2, "%llu\t %g\t %g \t %g\t %g %d\n", ihalo,(double)Radius,expected_R[0], expected_R[1],expected_R[2], block);
			    }
			}
		    }

		  /* calculate beta */
		  firstslope = (HaloTable[ihalo].Mvir - HaloTable[jhalo].Mvir) / (HaloTable[ihalo].Mvir + HaloTable[jhalo].Mvir);
		  firstslope /= (snapTime[iSnap] - snapTime[jSnap]) / (snapTime[iSnap] + snapTime[jSnap]);
		  firstslope = atan(firstslope)/asin(1.);
	
		  /* put it histogram */
		  block = (int) ((firstslope+1.)/binsize);
		  //printf("Block %d\n",block);
		  if(block == nStep) 
		    block = nStep-1;
		  /* discard when beta = 0 */
		  if(HaloTable[ihalo].Mvir < maxmass && HaloTable[ihalo].Mvir > minmass && HaloTable[jhalo].Mvir < maxmass && HaloTable[jhalo].Mvir > minmass)
		    {
		      if(SubTree[ihalo] == NULLPOINT && SubTree[jhalo] == NULLPOINT)
			{
			  MeanBeta[block]++;
			  if(block < 3 )
			    fprintf(fp3, "%llu\t%llu\n",i,ihalo);
			}
		    }

		  block = (int) ( (firstslope+1.) /binsize);

		  if(HaloTable[ihalo].Mvir < maxmass && HaloTable[ihalo].Mvir > minmass)
		    {
		      if(SubTree[ihalo] == NULLPOINT && SubTree[jhalo] == NULLPOINT)
			{
			  MeanBeta_h[block]++;
			}
		    }
		  /* calculate alpha */
		  firstslope =  (HaloTable[ihalo].Mvir - HaloTable[jhalo].Mvir)/(HaloTable[ihalo].Mvir + HaloTable[jhalo].Mvir);
		  /* put in histogram */
		  block = (int) ((firstslope+1.)/binsize);
		  //printf("Block %d\n",block);
		  if(block == nStep) 
		    block = nStep-1;
		  /* discard alpha = 0 */

		  if(HaloTable[ihalo].Mvir < maxmass && HaloTable[ihalo].Mvir > minmass)
		    Mean[block]++;

		  if(HaloTable[ihalo].Mvir < massboundery)
		    Mean_l[block]++;
		  else
		    Mean_h[block]++;
	  
		  //if(jhalo == 843642)
		  /* if jhalo has at least one progenitor */
		  if(output.progs[jhalo].nProgs > 0 )
		    {
		      /* set khalo = main progenitor of jhalo */
		      khalo = output.progs[jhalo].progID[0];
		      kSnap = HaloTable[khalo].SnapID;
		      if(HaloTable[khalo].Mvir > 1e-30)
			{
			  //printf("%llu %llu %llu\n",ihalo,jhalo,khalo);
			  //printf("%f %f %f\n",HaloTable[ihalo].Mvir,HaloTable[jhalo].Mvir,HaloTable[khalo].Mvir);

			  /* calculate beta_{n-1,n+1} */
			  firstslope = atan((HaloTable[ihalo].Mvir - HaloTable[khalo].Mvir) / (HaloTable[ihalo].Mvir + HaloTable[khalo].Mvir) /(snapTime[iSnap] - snapTime[kSnap]) * (snapTime[iSnap] + snapTime[kSnap]) ); 
			  /* calculate beta_{n-1,n} */
			  secondslope = atan((HaloTable[jhalo].Mvir - HaloTable[khalo].Mvir)/ (HaloTable[jhalo].Mvir + HaloTable[khalo].Mvir) /(snapTime[jSnap] - snapTime[kSnap]) * (snapTime[jSnap] + snapTime[kSnap]));

			  /* beta_{n-1,n} - beta_{n-1,n+1} */
			  corr_self = (secondslope-firstslope)/asin(1.)/2.;
	      

			  block = (int) ((corr_self + 1.)/binsize);
			  //printf("Block %d\n",block);
			  if(block == nStep) 
			    block = nStep-1;


			  if(HaloTable[ihalo].Mvir < maxmass && HaloTable[ihalo].Mvir > minmass && HaloTable[jhalo].Mvir < maxmass && HaloTable[jhalo].Mvir > minmass && HaloTable[khalo].Mvir < maxmass && HaloTable[khalo].Mvir > minmass )
			    {
			      if(SubTree[ihalo] == NULLPOINT && SubTree[jhalo] == NULLPOINT && SubTree[khalo] == NULLPOINT)
				CorrBeta[block]++;
			      
			    }
			  savedblock = block;
			  saveddouble = corr_self;

			  /* calculate alpha_{n-1,n+1} */
			  firstslope = atan((HaloTable[ihalo].Mvir - HaloTable[khalo].Mvir) / (lastmass) /(snapTime[iSnap] - snapTime[kSnap]) ); 
			  /* calculate alpha_{n-1,n} */
			  secondslope = atan((HaloTable[jhalo].Mvir - HaloTable[khalo].Mvir)/ (lastmass) /(snapTime[jSnap] - snapTime[kSnap]) );
			  /* alpha_{n-1,n} - alpha_{n-1,n+1} */
			  corr_self = (secondslope-firstslope)/asin(1.)/2.;
	      
			  block = (int) ((corr_self + 1.)/binsize);
			  //printf("Block %d\n",block);
			  if(block == nStep) 
			    block = nStep-1;
			  //printf("%f\t%d\n",corr_self,block);
			  
			  Corr[block]++;

			  if(HaloTable[jhalo].Mvir < massboundery)
			    Corr_l[block]++;
			  else
			    Corr_h[block]++;

			  R_k[0] = HaloTable[khalo].Xc;
			  R_k[1] = HaloTable[khalo].Yc;
			  R_k[2] = HaloTable[khalo].Zc;
			  
			  for(k=0;k<3;k++)
			    {
			      R1[k] = R_i[k] - R_j[k];
			      R2[k] = R_j[k] - R_k[k];
			      if(R1[k] > Boxsize/2.)
				R1[k] -= Boxsize;
			      else if(R1[k] < -1.*Boxsize/2)
				R1[k] += Boxsize;
			      if(R2[k] > Boxsize/2.)
				R2[k] -= Boxsize;
			      else if(R2[k] < -1.*Boxsize/2)
				R2[k] += Boxsize;				 
			    }
			  abs_R1 = sqrt(R1[0]*R1[0]+R1[1]*R1[1]+R1[2]*R1[2]);
			  abs_R2 = sqrt(R2[0]*R2[0]+R2[1]*R2[1]+R2[2]*R2[2]);
			  for(k=0;k<3;k++)
			    {
			      R1[k] /= abs_R1;
			      R2[k] /= abs_R2;
			    }
			  firstslope = R1[0]*R2[0]+R1[1]*R2[1]+R1[2]*R2[2];
			  block = (int) ((firstslope + 1.)/binsize);
			  if(block == nStep) 
			    block = nStep-1;
			  if(SubTree[ihalo] == NULLPOINT && SubTree[jhalo] == NULLPOINT && SubTree[khalo] == NULLPOINT)
			    {
			      postFluct[block]++;
			    }
			}
		    }
		  /* set ihalo = jhalo  */
		  ihalo = jhalo;
		}
	      else
		{
		  output.progs[ihalo].nProgs = 0;
		}
	    }
	  /* put count (link depth) in histogram */
	  //printf("count %d\n",count);
	  if(lastmass > minmass && lastmass < maxmass)
	    {
	      linkdepth[count]++;
	      if(count == 0)
		fprintf(fp1, "%llu\n", i);
	      if(SubTree[i] == NULLPOINT)
		{
		  hostlinkdepth[count]++;
		}
	      else
		{
		  resilinkdepth[count]++;
		}
	    }
	  /* put count in high/low mass histograms */
	  if(lastmass > massboundery)
	    {
	      highmasslinkdepth[count]++;
	    }
	  else
	    {
	      lowmasslinkdepth[count]++;
	    }

	}
    }

  printf("finish counting\n");
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);

  /* save outputs */
  sprintf(meanOutput,"%s_%dMeanSlope.dat",inputFile,highlim_npart);
  sprintf(l_meanOutput,"%s_%dMeanSlopeLowmass.dat",inputFile,highlim_npart);
  sprintf(h_meanOutput,"%s_%dMeanSlopeHighmass.dat",inputFile,highlim_npart);

  sprintf(meanbetaOutput,"%s_%dMeanBeta.dat",inputFile,highlim_npart);
  sprintf(l_meanbetaOutput,"%s_%dMeanBetaLowmass.dat",inputFile,highlim_npart);
  sprintf(h_meanbetaOutput,"%s_%dMeanBetaHighmass.dat",inputFile,highlim_npart);

  sprintf(corrOutput,"%s_%dSlopeCorr.dat",inputFile,highlim_npart);
  sprintf(l_corrOutput,"%s_%dSlopeCorrLowmass.dat",inputFile,highlim_npart);
  sprintf(h_corrOutput,"%s_%dSlopeCorrHighmass.dat",inputFile,highlim_npart);

  sprintf(corrbetaOutput,"%s_%dBetaCorr.dat",inputFile,highlim_npart);
  sprintf(l_corrbetaOutput,"%s_%dBetaCorrLowmass.dat",inputFile,highlim_npart);
  sprintf(h_corrbetaOutput,"%s_%dBetaCorrHighmass.dat",inputFile,highlim_npart);

  sprintf(depthOutput,"%s_%dLinkDepth.dat",inputFile,highlim_npart);
  sprintf(mainbranchmassOutput,"%s_%dMainbranchMass.dat",inputFile,highlim_npart);

  sprintf(dispOutput,"%s_%dDispMain.dat",inputFile,highlim_npart);
  sprintf(l_dispOutput,"%s_%dDispMainLowmass.dat",inputFile,highlim_npart);
  sprintf(h_dispOutput,"%s_%dDispMainHighmass.dat",inputFile,highlim_npart);

  sprintf(resi_dispOutput,"%s_%dDispMainResi.dat",inputFile,highlim_npart);
  sprintf(host_dispOutput,"%s_%dDispMainHost.dat",inputFile,highlim_npart);
  sprintf(postFluctOutput,"%s_%dPosFluctHost.dat", inputFile,highlim_npart);


  sprintf(alphasubOutput,"%s_%dAlphaSub.dat",inputFile,highlim_npart);
  sprintf(betasubOutput,"%s_%dBetaSub.dat",inputFile,highlim_npart);

  sprintf(lowmassLinkdepth,"%s_%dLinkDepthLowmass.dat",inputFile,highlim_npart);
  sprintf(highmassLinkdepth,"%s_%dLinkDepthHighmass.dat",inputFile,highlim_npart);  
  sprintf(host_Linkdepth,"%s_%dLinkDepthHost.dat",inputFile,highlim_npart); 
  sprintf(resi_Linkdepth,"%s_%dLinkDepthResi.dat",inputFile,highlim_npart);

  printf("finish make new file names\n");
  fp1 = fopen(meanOutput,"w+");
  fp2 = fopen(meanbetaOutput,"w+");
  fp3 = fopen(corrOutput, "w+");
  fp4 = fopen(corrbetaOutput, "w+");
  fp5 = fopen(dispOutput,"w+");

  printf("finish open files\n");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,Mean[i]);
      fprintf(fp2,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,MeanBeta[i]);
      fprintf(fp3,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,Corr[i]);
      fprintf(fp4,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,CorrBeta[i]);
      fprintf(fp5,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,dispCorr[i]);
    }
  printf("finish writing files\n");
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5);
  
  printf("finish closing files\n");

  free(Mean);
  //printf("mean\n");
  free(MeanBeta);
  //printf("mean beta\n");
  free(Corr);
  //printf("corr\n");
  free(CorrBeta);
  //printf("corr beta\n");
  free(dispCorr);
  //printf("distcorr\n");

  //printf("finish free parameters\n");
  //printf("finish lot 1\n");
  fp1 = fopen(l_meanOutput,"w+");
  fp2 = fopen(h_meanOutput,"w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,Mean_l[i]);
      fprintf(fp2,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,Mean_h[i]);
    }
  fclose(fp1);
  fclose(fp2);

  free(Mean_l);
  free(Mean_h);
  //printf("finish lot 2\n");
  fp1 = fopen(l_meanbetaOutput,"w+");
  fp2 = fopen(h_meanbetaOutput,"w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,MeanBeta_l[i]);
      fprintf(fp2,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,MeanBeta_h[i]);
    }
  fclose(fp1);
  fclose(fp2);

  free(MeanBeta_l);
  free(MeanBeta_h);
  //printf("finish lot 3\n");
  fp1 = fopen(l_corrOutput,"w+");
  fp2 = fopen(h_corrOutput,"w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,Corr_l[i]);
      fprintf(fp2,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,Corr_h[i]);
    }
  fclose(fp1);
  fclose(fp2);

  free(Corr_l);
  free(Corr_h);
  //printf("finish lot 4\n");
  fp1 = fopen(l_corrbetaOutput,"w+");
  fp2 = fopen(h_corrbetaOutput,"w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,CorrBeta_l[i]);
      fprintf(fp2,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,CorrBeta_h[i]);
    }
  fclose(fp1);
  fclose(fp2);

  free(CorrBeta_l);
  free(CorrBeta_h);
  //printf("finish lot 5\n");
  fp1 = fopen(l_dispOutput,"w+");
  fp2 = fopen(h_dispOutput,"w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,dispCorr_l[i]);
      fprintf(fp2,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,dispCorr_h[i]);
    }
  fclose(fp1);
  fclose(fp2);

  free(dispCorr_l);
  free(dispCorr_h);
  //printf("finish lot 6\n");
  fp1 = fopen(resi_dispOutput,"w+");
  fp2 = fopen(host_dispOutput,"w+");
  fp3 = fopen(postFluctOutput,"w+");
  for(i=0;i<nStep;i++)
    {
      fprintf(fp1,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,dispCorr_resi[i]);
      fprintf(fp2,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,dispCorr_host[i]);
      fprintf(fp3,"%f\t%llu\n",-1.+((float)(i)+.5)*binsize,postFluct[i]);
    }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  free(dispCorr_resi);
  free(dispCorr_host);
  free(postFluct);
  //printf("finish lot 7\n");
  fp1 = fopen(depthOutput,"w+");
  fp2 = fopen(mainbranchmassOutput,"w+");
  fp3 = fopen(highmassLinkdepth,"w+");
  fp4 = fopen(lowmassLinkdepth,"w+");
  fp5 = fopen(host_Linkdepth,"w+");
  fp6 = fopen(resi_Linkdepth,"w+");
  for(i=0;i<NSNAPS;i++)
    {
      fprintf(fp1,"%llu\t%llu\n",i,linkdepth[i]);
      fprintf(fp4,"%llu\t%llu\n",i,lowmasslinkdepth[i]);
      fprintf(fp3,"%llu\t%llu\n",i,highmasslinkdepth[i]);
      fprintf(fp5,"%llu\t%llu\n",i,hostlinkdepth[i]);
      fprintf(fp6,"%llu\t%llu\n",i,resilinkdepth[i]);
      fprintf(fp2,"%g\t%g\n",snapTime[i],meanmainbranch[i]/ (float)TotMainbranch[i]);
    }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5);
  fclose(fp6);
  free(hostlinkdepth);
  free(linkdepth);
  free(lowmasslinkdepth);
  free(highmasslinkdepth);
  free(meanmainbranch);
  free(TotMainbranch);
  //printf("finish lot 8\n");
  printf("finish subroutine\n");
}

