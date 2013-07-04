#include "allvars.h"

MyIDtype findtoplevel(MyIDtype hid)
{
  MyIDtype hostid;
  hostid = SubTree[hid];
  if(hostid < NULLPOINT)
    {
      hostid = findtoplevel(hostid);
    }
  else
    return hid;
}

void printoutfullAHF()
{
  MyIDtype ihalo;
  char header_out[MAXSTRING];
  sprintf(header_out, 
	  "#ID(1)     hostHalo(2)     numSubStruct(3) Mvir(4) npart(5)        Xc(6)   Yc(7)   Zc(8)   VXc(9)  VYc(10) VZc(11) Rvir(12)        Rmax(13)        r2(14)  mbp_offset(15)  com_offset(16)  Vmax(17)        v_esc(18)       sigV(19)        lambda(20)      lambdaE(21)     Lx(22)  Ly(23)  Lz(24)  b(25)   c(26)   Eax(27) Eay(28) Eaz(29) Ebx(30) Eby(31) Ebz(32) Ecx(33) Ecy(34) Ecz(35) ovdens(36)      nbins(37)       fMhires(38)     Ekin(39)        Epot(40)        SurfP(41)       Phi0(42)        cNFW(43)");
  //printf("%s\n",header_out);
  printf("%s\n", header_out);
  for(ihalo=0;ihalo<TotNhalos;ihalo++)
    {
  printf("%llu\t%llu\t%lu\t%.8g\t%lu\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\n",
	  HaloTable[ihalo].AHFID,
	  HaloTable[ihalo].hostHalo,
	  HaloTable[ihalo].numSubStruct,
	  HaloTable[ihalo].Mvir,
	  HaloTable[ihalo].npart,
	  HaloTable[ihalo].Xc,
	  HaloTable[ihalo].Yc,
	  HaloTable[ihalo].Zc,
	  HaloTable[ihalo].VXc,
	  HaloTable[ihalo].VYc,
	  HaloTable[ihalo].VZc,
	  HaloTable[ihalo].Rvir,
	  HaloTable[ihalo].Rmax,
	  HaloTable[ihalo].r2,
	  HaloTable[ihalo].mbp_offset,
	  HaloTable[ihalo].com_offset,
	  HaloTable[ihalo].Vmax,
	  HaloTable[ihalo].v_esc,
	  HaloTable[ihalo].sigV,
	  HaloTable[ihalo].lambda,
	  HaloTable[ihalo].lambdaE,
	  HaloTable[ihalo].Lx,
	  HaloTable[ihalo].Ly,
	  HaloTable[ihalo].Lz,
	  HaloTable[ihalo].b,
	  HaloTable[ihalo].c,
	  HaloTable[ihalo].Eax,
	  HaloTable[ihalo].Eay,
	  HaloTable[ihalo].Eaz,
	  HaloTable[ihalo].Ebx,
	  HaloTable[ihalo].Eby,
	  HaloTable[ihalo].Ebz,
	  HaloTable[ihalo].Ecx,
	  HaloTable[ihalo].Ecy,
	  HaloTable[ihalo].Ecz,
	  HaloTable[ihalo].ovdens,
	  HaloTable[ihalo].nbins,
	  HaloTable[ihalo].fMhires,
	  HaloTable[ihalo].Ekin,
	  HaloTable[ihalo].Epot,
	  HaloTable[ihalo].SurfP,
	  HaloTable[ihalo].Phi0,
	  HaloTable[ihalo].cNFW 
	  );
    }
}



void makesubfindout()
{
  FILE *fp1, *fp2, *fp3, *fp4;
  struct subfind_data data;
  MyIDtype *id_to_sub;
  MyIDtype ngroups,nids,nsubs,firstsub, lastsub;
  MyIDtype hostid,groupid,subid,*subhalolen;
  MyIDtype iF,j,i,k,pid,hid,ahfid,countsub,countgroup,countpid,startj,stopj,startsub,startpid;
  char buf[MAXSTRING], buf2[MAXSTRING];

  int bufferint;
  float bufferfloat;
  long long bufferlong;
  MyIDtype **grouplist;
  MyIDtype *sub_to_group;
  printf("Start making subfind outputs\n");
  countsub = 0;
  sprintf(buf, "%s/SUBFIND/Mapback_order.txt", FolderName);      
  fp3 = fopen(buf,"w+");
  sprintf(buf, "%s/SUBFIND/Mapback_hash.txt", FolderName);      
  fp4 = fopen(buf,"w+");
  for(iF=FIRSTSNAP;iF<=LASTSNAP;iF++)
    {
            
      sprintf(buf2, "mkdir -p %s/SUBFIND/groups_%03d/",FolderName, (int)iF);
      system(buf2);
      sprintf(buf, "%s/SUBFIND/groups_%03d/subhalo_tab_%03d.0", FolderName, (int)iF, (int)iF);
      fp1 = fopen(buf,"wb");
      sprintf(buf, "%s/SUBFIND/groups_%03d/subhalo_ids_%03d.0", FolderName, (int)iF, (int)iF);      
      fp2 = fopen(buf,"wb");

      id_to_sub = malloc((box_particle+1)*sizeof(MyIDtype));
      for(i=0;i<=box_particle;i++)
	id_to_sub[i] = NULLPOINT;
      if(SnapNhalos[iF] == 0)
	{
	  data.TotNgroups = 0;
	  data.TotNids = 0;
	  data.TotNsubhalos = 0;
	  

	  fwrite(&(data.TotNgroups), sizeof(int), 1, fp1);
	  fwrite(&(data.TotNgroups), sizeof(int), 1, fp1);
	  fwrite(&(data.TotNids), sizeof(int), 1, fp1);
	  fwrite(&(data.TotNids), sizeof(long long), 1, fp1);
	  bufferint = 1;
	  fwrite(&(bufferint), sizeof(int), 1, fp1);
	  fwrite(&(data.TotNsubhalos), sizeof(int), 1, fp1);
	  fwrite(&(data.TotNsubhalos), sizeof(int), 1, fp1);

	  fwrite(&(data.TotNgroups), sizeof(int), 1, fp2);
	  fwrite(&(data.TotNgroups), sizeof(int), 1, fp2);
	  fwrite(&(data.TotNids), sizeof(int), 1, fp2);
	  fwrite(&(data.TotNids), sizeof(long long), 1, fp2);
	  bufferint = 1;
	  fwrite(&(bufferint), sizeof(int), 1, fp2);
	  bufferint = 0;
	  // set offset = 0
	  fwrite(&(bufferint), sizeof(int), 1, fp2);


	}
      else
	{
	  firstsub = countsub;
	  lastsub = firstsub + SnapNhalos[iF] - 1;
	  printf("%llu first:%llu last:%llu\n",iF,firstsub,lastsub);
	  ngroups = 0;
	  nsubs = SnapNhalos[iF];
	  printf("allocate id_to_sub\n");
     
	  printf("loop for ngroups and id_to_sub\n");

	  for(j=firstsub; j<= lastsub; j++)
	    {
	      //printf("j = %d\n",j);
	      hid = j-firstsub;
	      if(SubTree[j] == NULLPOINT)
		{
		  ngroups ++;
		}	    
	      for(i=0;i<HaloTable[j].npart; i++)
		{
		  pid = HaloTable[j].Particles[i].ParticleID;
		  //printf("pid = %llu\n",pid);
	
		  id_to_sub[pid] = hid;
		}
	    }
	  nids = 0;
	  printf("loop to count nids\n");
	  for(i=0;i<=box_particle;i++)
	    {
	      if(id_to_sub[i] < NULLPOINT)
		nids++;
	    }
	  data.TotNgroups = (int)ngroups;
	  data.TotNids = (long long)nids;
	  data.TotNsubhalos =(int) nsubs;

	  fwrite(&(data.TotNgroups), sizeof(int), 1, fp1);
	  fwrite(&(data.TotNgroups), sizeof(int), 1, fp1);
	  fwrite(&(data.TotNids), sizeof(int), 1, fp1);
	  fwrite(&(data.TotNids), sizeof(long long), 1, fp1);
	  bufferint = 1;
	  fwrite(&(bufferint), sizeof(int), 1, fp1);
	  fwrite(&(data.TotNsubhalos), sizeof(int), 1, fp1);
	  fwrite(&(data.TotNsubhalos), sizeof(int), 1, fp1);

	  fwrite(&(data.TotNgroups), sizeof(int), 1, fp2);
	  fwrite(&(data.TotNgroups), sizeof(int), 1, fp2);
	  fwrite(&(data.TotNids), sizeof(int), 1, fp2);
	  fwrite(&(data.TotNids), sizeof(long long), 1, fp2);
	  bufferint = 1;
	  fwrite(&(bufferint), sizeof(int), 1, fp2);
	  bufferint = 0;
	  // set offset = 0
	  fwrite(&(bufferint), sizeof(int), 1, fp2);



	  printf("allocate memory\n");
	  data.IdList = malloc(sizeof(unsigned int) *data.TotNids);
	  data.SubLen = calloc(data.TotNsubhalos,sizeof(int) );
	  data.SubOffset = malloc(sizeof(unsigned int) *data.TotNsubhalos);

	  data.SubParentHalo = malloc(sizeof(int) *data.TotNsubhalos);
	  data.SubhaloMass = malloc(sizeof(float) *data.TotNsubhalos);

	  data.SubhaloPos = malloc(3 * sizeof(float ) *data.TotNsubhalos);
	  data.SubhaloVel = malloc(3 * sizeof(float ) *data.TotNsubhalos);
	  data.SubhaloCM = malloc(3 * sizeof(float ) *data.TotNsubhalos);
	  data.SubhaloSpin = malloc(3 * sizeof(float ) *data.TotNsubhalos);


	  data.SubhaloVelDisp = malloc(sizeof(float) *data.TotNsubhalos);
	  data.SubhaloVmax = malloc(sizeof(float) *data.TotNsubhalos);
	  data.SubhaloVmaxRad = malloc(sizeof(float) *data.TotNsubhalos);
	  data.SubhaloHalfMass = malloc(sizeof(float) *data.TotNsubhalos);

	  data.SubhaloMostBoundID = malloc(sizeof(unsigned int) *data.TotNsubhalos);
	  data.SubhaloGrNr = malloc(sizeof(int) *data.TotNsubhalos);

	  data.GroupNsubs = calloc(data.TotNgroups,sizeof(int));

	  data.GroupMass = malloc(sizeof(float) *data.TotNgroups);
	  data.GroupPos = malloc(3 * sizeof(float) *data.TotNgroups);

	  data.GroupOffset = malloc(data.TotNgroups*sizeof(unsigned int));
	  data.GroupLen = calloc(data.TotNgroups,sizeof(int));
	  data.Group_M_Mean200 = malloc(sizeof(float) *data.TotNgroups);
	  data.Group_R_Mean200 = malloc(sizeof(float) *data.TotNgroups);
	  data.Group_M_Crit200 = malloc(sizeof(float) *data.TotNgroups);
	  data.Group_R_Crit200 = malloc(sizeof(float) *data.TotNgroups);
	  data.GroupContaminationCount = malloc(sizeof(int) *data.TotNgroups);
	  data.GroupContaminationMass = malloc(sizeof(float) *data.TotNgroups);

	  data.Group_M_TopHat200 = malloc(sizeof(float) *data.TotNgroups);
	  data.Group_R_TopHat200 = malloc(sizeof(float) *data.TotNgroups);

	  data.GroupFirstSub = malloc(sizeof(int) *data.TotNgroups);
	  grouplist = malloc(data.TotNgroups*sizeof(MyIDtype));

	  sub_to_group = malloc(data.TotNsubhalos*sizeof(MyIDtype));

	  for(j=0;j<data.TotNgroups;j++)
	    {
	      grouplist[j] = malloc(0*sizeof(MyIDtype));
	    }

	  for(j=firstsub; j<= lastsub; j++)
	    {
	      sub_to_group[j-firstsub] = findtoplevel(j)-firstsub;
	      
	      //printf("%llu => %llu\n",j-firstsub,sub_to_group[j-firstsub]);
	    }
	  countgroup = 0;
	  for(j=0; j< SnapNhalos[iF]; j++)
	    {
	      if(sub_to_group[j] == j)
	  	{
	  	  sub_to_group[j] = countgroup;
	  	  countgroup++;
	  	}
	      else
	  	{
	  	  sub_to_group[j] = sub_to_group[sub_to_group[j]];
	  	}
	      //printf("%llu => %llu\n",j,sub_to_group[j]);
	    }
	  if(countgroup != data.TotNgroups)
	    {
	      printf("Subgroup doesn't match\n");
	      exit(0);
	    }
	  
	  /* get GroupNsubs and lists */
	  printf("get GroupNsubs and lists\n");
	  for(j=0; j< data.TotNsubhalos; j++)
	    {
	      groupid = sub_to_group[j];
	      data.GroupNsubs[groupid]++;
	      grouplist[groupid] = realloc(grouplist[groupid],sizeof(MyIDtype)*data.GroupNsubs[groupid]);
	      //printf("groupid:%llu N:%d\n",groupid,data.GroupNsubs[groupid]);
	      grouplist[groupid][data.GroupNsubs[groupid]-1] = j;
	    }
	  
	  free(sub_to_group);

	  /* get SubLen & GroupLen & IDlist*/
	  printf("Get SubLen\n");
	  
	  subhalolen = calloc(data.TotNsubhalos,sizeof(MyIDtype));
	  for(i=0;i<=box_particle;i++)
	    {
	      subid = id_to_sub[i];
	      if(subid < NULLPOINT)
	  	subhalolen[subid]++;
	    }

	  printf("Get GroupLen\n");
	  startsub = 0;
	  startpid = 0;
	  for(j=0;j<data.TotNgroups;j++)
	    {
	      data.GroupLen[j] = 0;
	      for(i=0;i<data.GroupNsubs[j];i++)
	      	{
	      	  subid = grouplist[j][i];
	      	  ahfid = subid + firstsub;
	      	  //printf("AHF id:%llu\n",ahfid);
	      	  data.SubLen[startsub] = (int) subhalolen[subid];
	      	  data.GroupLen[j] += (int) subhalolen[subid];
		  data.SubParentHalo[startsub] = SubTree[ahfid];

		  if(data.SubParentHalo[startsub] == NULLPOINT)
		    {
		      data.SubParentHalo[startsub] = 0;
		    }
		  else
		    {
		      hostid = data.SubParentHalo[startsub] - firstsub;
		      for(k=0;k<data.GroupNsubs[j];k++)
			{
			  if(hostid == grouplist[j][k])
			    {
			      data.SubParentHalo[startsub] = (int) k;
			    }
			}
		    }
#ifdef STORESUBIDINMOSTBOUNDID
		  data.SubhaloMostBoundID[j] = (unsigned int) ahfid;
		  fprintf(fp4,"%llu\t%llu\n",ahfid, HaloTable[ahfid].ID);
#endif
		  fprintf(fp3,"%llu\t%llu\n",MAXHALOPERSNAP*iF + startsub +1, HaloTable[ahfid].ID);

		  //fprintf(fp3,"%llu\t%llu\n",MAXHALOPERSNAP*iF + startsub +1, HaloTable[ahfid].ID);
		  data.SubhaloGrNr[startsub] = (int) j;
	
		  data.SubhaloMass[startsub] = (float) (data.SubLen[startsub]*p_mass/GagetUnit2Msun);
		  
		  data.SubhaloPos[3*startsub] = (float) (HaloTable[ahfid].Xc/Mpc2kpc);
		  data.SubhaloPos[3*startsub+1] = (float) (HaloTable[ahfid].Yc/Mpc2kpc);
		  data.SubhaloPos[3*startsub+2] = (float) (HaloTable[ahfid].Zc/Mpc2kpc);

		  data.SubhaloVel[3*startsub] = (float) HaloTable[ahfid].VXc;
		  data.SubhaloVel[3*startsub+1] = (float) HaloTable[ahfid].VYc;
		  data.SubhaloVel[3*startsub+2] = (float) HaloTable[ahfid].VZc;

		  /* skip SubhaloCM */
		  data.SubhaloCM[3*startsub] = (float) default_float;
		  data.SubhaloCM[3*startsub+1] = (float) default_float;
		  data.SubhaloCM[3*startsub+2] = (float) default_float;

		  data.SubhaloSpin[3*startsub] = (float) HaloTable[ahfid].Lx;
		  data.SubhaloSpin[3*startsub+1] = (float) HaloTable[ahfid].Ly;
		  data.SubhaloSpin[3*startsub+2] = (float) HaloTable[ahfid].Lz;

		  data.SubhaloVelDisp[startsub] =  (float) HaloTable[ahfid].sigV;
		  data.SubhaloVmax[startsub] = (float) HaloTable[ahfid].Vmax;
		  data.SubhaloVmaxRad[startsub] = (float) HaloTable[ahfid].Rmax;

		  /* skip SubhaloHalfMass */
		  data.SubhaloHalfMass[startsub] = (float) default_float;

	      	  for(k=0;k<HaloTable[ahfid].npart;k++)
	      	    {
	      	      pid = HaloTable[ahfid].Particles[k].ParticleID;
	      	      if(id_to_sub[pid] == subid)
	      		{
	      		  data.IdList[startpid] = (unsigned int) pid;
	      		  startpid++;
	      		}
	      	    }
		  startsub++;
	      	}
	      //printf("Group:%llu %d pids\n",j,data.GroupLen[j]);
	    }
	  free(subhalolen);
	  /* Get GroupOffset */
	  data.GroupOffset[0] = 0;
	  data.GroupFirstSub[0] = 0;
	  for(j=1;j<data.TotNgroups;j++)
	    {
	      data.GroupOffset[j] = data.GroupOffset[j-1] + data.GroupLen[j-1];
	      data.GroupFirstSub[j] =  data.GroupFirstSub[j-1] + data.GroupNsubs[j-1];
	    }
	  /* Get Suboffset */
	  data.SubOffset[0] = 0;
	  for(j=1;j<data.TotNsubhalos;j++)
	    {
	      data.SubOffset[j] =  (int) (data.SubOffset[j-1]+data.SubLen[j-1]);
	      printf("%d : %lu / %lu => %lu / %lu\n",iF,j,data.TotNsubhalos,data.SubOffset[j],data.TotNids);
	    }
#ifndef STORESUBIDINMOSTBOUNDID
	  for(j=0;j<data.TotNsubhalos;j++)
	    {
	      data.SubhaloMostBoundID[j] = data.IdList[data.SubOffset[j]];
	    }
#endif
	  /* Get Main halo properties */
	  printf("Get halo's main properties\n");
	  for(j=0;j<data.TotNgroups;j++)
	    {
	      subid = grouplist[j][0] + firstsub;
	      data.GroupMass[j] = (float) (data.GroupLen[j]*p_mass/GagetUnit2Msun);
	      data.GroupPos[3*j] = (float) (HaloTable[subid].Xc/Mpc2kpc);
	      data.GroupPos[3*j+1] = (float) (HaloTable[subid].Yc/Mpc2kpc);
	      data.GroupPos[3*j+2] = (float) (HaloTable[subid].Zc/Mpc2kpc);
	      data.Group_M_Mean200[j] = (float) (HaloTable[subid].Mvir/GagetUnit2Msun);
	      data.Group_R_Mean200[j] = (float) (HaloTable[subid].Rvir);
	      /* skip  Group_M_Crit200 */
	      data.Group_M_Crit200[j] = (float) default_float;
	      /* skip  Group_R_Crit200 */
	      data.Group_R_Crit200[j] = (float) default_float;
	      /* skip  Group_M_TopHat200 */
	      data.Group_M_TopHat200[j] = (float) default_float;
	      /* skip  Group_R_TopHat200 */
	      data.Group_R_TopHat200[j] = (float) default_float;
	      /*  skip  GroupContaminationCount  */
	      data.GroupContaminationCount[j] = (int) default_int;
	      /* skip  GroupContaminationMass */
	      data.GroupContaminationMass[j] = (float) default_int;
	    }

		 printf("Start writing\n");	
	  /* Start writing */
	  fwrite(&(data.GroupLen[0]), sizeof(int), data.TotNgroups, fp1);
	  fwrite(&(data.GroupOffset[0]), sizeof(unsigned int), data.TotNgroups, fp1);
	  fwrite(&(data.GroupMass[0]), sizeof(float), data.TotNgroups, fp1);
	  fwrite(&(data.GroupPos[0]), 3*sizeof(float), data.TotNgroups, fp1);

	  fwrite(&(data.Group_M_Mean200[0]), sizeof(float), data.TotNgroups, fp1);
	  fwrite(&(data.Group_R_Mean200[0]), sizeof(float), data.TotNgroups, fp1);
	  fwrite(&(data.Group_M_Crit200[0]), sizeof(float), data.TotNgroups, fp1);
	  fwrite(&(data.Group_R_Crit200[0]), sizeof(float), data.TotNgroups, fp1);
	  fwrite(&(data.Group_M_TopHat200[0]), sizeof(float), data.TotNgroups, fp1);
	  fwrite(&(data.Group_R_TopHat200[0]), sizeof(float), data.TotNgroups, fp1); 

#ifdef SO_VEL_DISPERSIONS
/* 	  fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_VelDisp_Mean200 *\/ */
/* 	  fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_VelDisp_Crit200 *\/ */
/* 	  fseek(fd, sizeof(float) * ngroups, SEEK_CUR);	/\* skip  Group_VelDisp_TopHat200 *\/ */
#endif
	  fwrite(&(data.GroupContaminationCount[0]), sizeof(int), data.TotNgroups, fp1);
	  fwrite(&(data.GroupContaminationMass[0]), sizeof(float), data.TotNgroups, fp1);
	  fwrite(&(data.GroupNsubs[0]), sizeof(int), data.TotNgroups, fp1); 
	  fwrite(&(data.GroupFirstSub[0]), sizeof(int), data.TotNgroups, fp1);
	  fwrite(&(data.SubLen[0]), sizeof(int), data.TotNsubhalos, fp1);

	  fwrite(&(data.SubOffset[0]), sizeof(int), data.TotNsubhalos, fp1);


	  fwrite(&(data.SubParentHalo[0]), sizeof(int),data.TotNsubhalos , fp1);
	  fwrite(&(data.SubhaloMass[0]), sizeof(float),data.TotNsubhalos , fp1);

	  fwrite(&(data.SubhaloPos[0]), 3 * sizeof(float),data.TotNsubhalos , fp1);
	  fwrite(&(data.SubhaloVel[0]), 3 * sizeof(float),data.TotNsubhalos , fp1); 
	  fwrite(&(data.SubhaloCM[0]), 3 * sizeof(float),data.TotNsubhalos , fp1);
	  fwrite(&(data.SubhaloSpin[0]), 3 * sizeof(float),data.TotNsubhalos , fp1);

	  fwrite(&(data.SubhaloVelDisp[0]), sizeof(float),data.TotNsubhalos , fp1);
	  fwrite(&(data.SubhaloVmax[0]), sizeof(float),data.TotNsubhalos , fp1);
	  fwrite(&(data.SubhaloVmaxRad[0]), sizeof(float),data.TotNsubhalos , fp1);

	  fwrite(&(data.SubhaloHalfMass[0]), sizeof(float),data.TotNsubhalos , fp1);
	  fwrite(&(data.SubhaloMostBoundID[0]), sizeof(unsigned int),data.TotNsubhalos , fp1);
	  fwrite(&(data.SubhaloGrNr[0]), sizeof(int), data.TotNsubhalos, fp1);

	  /* write particle ids */
	  fwrite(&(data.IdList[0]), sizeof(unsigned int), data.TotNids, fp2);
	  free(data.IdList);
 
	  free(data.SubLen);
	  free(data.SubOffset);

	  free(data.SubParentHalo);
	  free(data.SubhaloMass);

	  free(data.SubhaloPos);
	  free(data.SubhaloVel);
	  free(data.SubhaloCM);
	  free(data.SubhaloSpin);


	  free(data.SubhaloVelDisp);
	  free(data.SubhaloVmax);
	  free(data.SubhaloVmaxRad);
	  free(data.SubhaloHalfMass);

	  free(data.SubhaloMostBoundID);
	  free(data.SubhaloGrNr);

	  free(data.GroupNsubs);
	  free(data.GroupMass);
	  free(data.GroupPos);
	  
	  free(data.GroupOffset);
	  free(data.GroupLen);
	  free(data.Group_M_Mean200);
	  free(data.Group_R_Mean200);
	  free(data.Group_M_Crit200);
	  free(data.Group_R_Crit200);
	  free(data.Group_M_TopHat200);
	  free(data.Group_R_TopHat200);

	  free(data.GroupContaminationCount);
	  free(data.GroupContaminationMass);
	  free(data.GroupFirstSub);

	  free(id_to_sub);

	  for(j=0;j<data.TotNgroups;j++)
	    {
	      free(grouplist[j]);
	    }
	  free(grouplist);
	}

      
      
      fclose(fp1);
      fclose(fp2);
      printf("Closing files\n");
      countsub += SnapNhalos[iF];
    }
  fclose(fp3);
  fclose(fp4);
  printf("Finish subfindoutput()\n");
  exit(0);
} 
void printoutprofile()
{
  FILE *fp,*bfp;
  char filename[MAXSTRING];
  MyIDtype i;
  sprintf(filename,"%s/%s.profile",TreeFolder,inputFile);
  fp = fopen(filename, "w+");
  sprintf(filename,"%s/%s.profile.bin",TreeFolder,inputFile);
  bfp = fopen(filename, "wb");
  for(i=0;i<TotNhalos;i++)
    {
      fprintf(fp,"%llu\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	      HaloTable[i].ID,
	      HaloTable[i].SnapID,
	      HaloTable[i].Mvir,
	      HaloTable[i].Rvir,
	      HaloTable[i].Xc,
	      HaloTable[i].Yc,
	      HaloTable[i].Zc,
	      HaloTable[i].VXc,
	      HaloTable[i].VYc,
	      HaloTable[i].VZc);
      fwrite (&(HaloTable[i].ID) , 1 , sizeof(MyIDtype) , bfp );
      fwrite (&(HaloTable[i].SnapID) , 1 , sizeof(unsigned int) , bfp );
      fwrite (&(HaloTable[i].Mvir) , 1 , sizeof(float) , bfp );
      fwrite (&(HaloTable[i].Rvir) , 1 , sizeof(float) , bfp );
      fwrite (&(HaloTable[i].Xc) , 1 , sizeof(float) , bfp );
      fwrite (&(HaloTable[i].Yc) , 1 , sizeof(float) , bfp );
      fwrite (&(HaloTable[i].Zc) , 1 , sizeof(float) , bfp );
      fwrite (&(HaloTable[i].VXc) , 1 , sizeof(float) , bfp );
      fwrite (&(HaloTable[i].VYc) , 1 , sizeof(float) , bfp );
      fwrite (&(HaloTable[i].VZc) , 1 , sizeof(float) , bfp );
    }
  fclose(fp);
  fclose(bfp);
}

void printoutparticles_binary()
{
  FILE *fp_npart,*fp_particle, *fp_uid;
  char filename[MAXSTRING];
  MyIDtype i,j;
  if(output.outputFormat > 1.129 && output.outputFormat < 1.131)
    {
      sprintf(filename,"%s/%s_particles.bin",FolderName,inputFile);
      fp_particle = fopen(filename, "wb");
      sprintf(filename,"%s/%s_npart.bin",FolderName,inputFile);
      fp_npart = fopen(filename, "wb");
      sprintf(filename,"%s/%s_uid.bin",FolderName,inputFile);
      fp_uid = fopen(filename, "wb");
    }
  else
    {
      sprintf(filename,"%s/all_particles.bin",FolderName);
      fp_particle = fopen(filename, "wb");
      sprintf(filename,"%s/all_npart.bin",FolderName);
      fp_npart = fopen(filename, "wb");
      sprintf(filename,"%s/all_uid.bin",FolderName);
      fp_uid = fopen(filename, "wb");
    }


  fwrite (&TotNhalos , 1 , sizeof(MyIDtype) , fp_uid);
  fwrite (&TotNhalos , 1 , sizeof(MyIDtype) , fp_npart);
  fwrite (&TotNhalos , 1 , sizeof(MyIDtype) , fp_particle);
  for(i=0;i<TotNhalos;i++)
    {
      fwrite (&(HaloTable[i].ID) , 1 , sizeof(MyIDtype) , fp_uid );
      fwrite (&(HaloTable[i].npart) , 1 , sizeof(MyIDtype) , fp_npart );
      for(j=0;j<HaloTable[i].npart;j++)
	{
	  fwrite (&(HaloTable[i].Particles[j].ParticleID) , 1 , sizeof(MyIDtype), fp_particle);
	}
    }
  fclose(fp_uid);
  fclose(fp_npart);
  fclose(fp_particle);
}


void printouttrees()
{
  FILE* fp;
  char filename[MAXSTRING];
  MyIDtype i,j;
  sprintf(filename,"%s/%s.trees",TreeFolder,inputFile);
  fp = fopen(filename, "w+");
  for(i=0;i<TotNhalos;i++)
    {
      fprintf(fp, "%llu\t", i);
      for(j=0;j<output.progs[i].nProgs; j++)
	{
	  fprintf(fp,"%llu\t",output.progs[i].progID[j]);
	}
      fprintf(fp,"\n");
    }
  fclose(fp);
}
