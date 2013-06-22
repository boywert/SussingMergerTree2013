#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <Python.h>
#include <numpy/arrayobject.h>

//global vars -- not very good, but I do not want to fix it at the moment
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

//python main function
static PyObject *
gadget(PyObject *self, PyObject *args)
{
  const char *file;
  int i,j,k,pos,nfiles,tempint,dummy,ntot_withmasses;
  int NumPart,Ngas,t,n,off,pc,pc_new,pc_sph;
  npy_intp dim[1] = {1};
  npy_intp ndim = 1;
  PyArrayObject *array;
  float *temp,tempfloat;
  struct io_header header;
  FILE *fd;
  char buf[500];
  
  if (!PyArg_ParseTuple(args, "si", &file, &nfiles))  
    return NULL;

  
  //read from file(s) in a loop
  //***keep track of global particle index through pc and local file index through pc_new***
  pc = 0;
  pc_new = 0;
  for(i=0;i<nfiles;i++)
    {
      pc = pc_new;
      if(nfiles > 1)  //multiple file names are different
	sprintf(buf,"%s.%d",file,i);
      else
	sprintf(buf,"%s",file);

      if(!(fd=fopen(buf,"r")))  //error check and complain
	{
	  printf("gadgetIO module can't open file: '%s'\n",buf);
	  fflush(stdout);
	  return NULL;
	}

      printf("reading: '%s'\n",buf); 
      fflush(stdout);
      
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header, sizeof(header), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);
      
      //allocate data array and get total number of particles to start reading
      if(i == 0)
	{
	  if(nfiles == 1)
	    {
	      for(k=0,NumPart=0;k<5;k++)
		NumPart += header.npart[k];
	      Ngas = header.npart[0];
	    }
	  else
	    {
	      for(k=0,NumPart=0;k<5;k++)
		NumPart += header.npartTotal[k];
	      Ngas = header.npartTotal[0];
	    }
	  
	  //make numpy array and copy info into it
	  dim[0] = NumPart*3;
	  array = (PyArrayObject *) PyArray_SimpleNew(ndim, dim, NPY_FLOAT);
	}
      
      for(k=0,ntot_withmasses=0;k<5;k++)
	if(header.mass[k] == 0)
	  ntot_withmasses += header.npart[k];

      fread(&dummy, sizeof(dummy), 1, fd);
      for(k=0;k<6;k++)
	{ 
	  printf("Read particle %d : N = %d\n",k,header.npart[k]);
	  fread(PyArray_GETPTR1(array,pc_new),sizeof(float),header.npart[k]*3,fd);
	  pc_new += header.npart[k]*3;
	}
      fread(&dummy, sizeof(dummy), 1, fd);
      fclose(fd);
    }
  
  array = (PyArrayObject *) PyArray_Return(array);
  
  return (PyObject *) array;
}

PyMethodDef methods[] = {
  {"readone", gadget, METH_VARARGS, "really really really crappy code to read in gadget positions and velocities...enjoy! -Matthew Becker UofC 2009 \n x = gadgetPyIO.readone(filename,nfiles,0) "},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC 
initgadgetPyIO()
{
  (void) Py_InitModule("gadgetPyIO", methods);   
  import_array();
}

/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 int unit_conversion(void)
 {
 double GRAVITY, BOLTZMANN, PROTONMASS;
 double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
 double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;  
 double G, Xh, HubbleParam;
 
 int i;
 double MeanWeight, u, gamma;
 
 //physical constants in cgs units
 GRAVITY   = 6.672e-8;
 BOLTZMANN = 1.3806e-16;
 PROTONMASS = 1.6726e-24;
 
 //internal unit system of the code
 UnitLength_in_cm= 3.085678e21;   //code length unit in cm/h
 UnitMass_in_g= 1.989e43;         //code mass unit in g/h
 UnitVelocity_in_cm_per_s= 1.0e5;
 
 UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
 UnitDensity_in_cgs= UnitMass_in_g/ pow(UnitLength_in_cm,3);
 UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ pow(UnitTime_in_s,2);
 UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
 
 G=GRAVITY/ pow(UnitLength_in_cm,3) * UnitMass_in_g * pow(UnitTime_in_s,2);
 
 
 Xh= 0.76;  //mass fraction of hydrogen
 HubbleParam= 0.65;
 
 for(i=1; i<=NumPart; i++)
 {
 if(P[i].Type==0)  //gas particle
 {
 MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].Ne) * PROTONMASS;
 
 //convert internal energy to cgs units
 
 u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;
 
 gamma= 5.0/3;
 
 //get temperature in Kelvin
 
 P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
 }
 }
 }
*/
