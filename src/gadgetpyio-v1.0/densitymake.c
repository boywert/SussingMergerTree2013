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
  PyObject *array;
  float *data;
  float *temp,tempfloat;
  struct io_header header;
  FILE *fd;
  float *dataout;
  char buf[500];
  double x0,y0,z0,xmin,xmax,ymin,ymax,cellsize,rho,boxsize,radius,length;
  double x_axis[3],y_axis[3],z_axis[3],rs[0],x,y,double_buf;
  int tot_x, tot_y, x_block,y_block;
  double **count;
  long countp,countk;
  size_t size,il;

  if (!PyArg_ParseTuple(args, "siddddddddddddddddddii", 
			&file, &nfiles,
			&x0,&y0,&z0,
			&(x_axis[0]),&(x_axis[1]),&(x_axis[2]),
			&(y_axis[0]),&(y_axis[1]),&(y_axis[2]),
			&(z_axis[0]),&(z_axis[1]),&(z_axis[2]),
			&xmin,&xmax,
			&ymin,&ymax,
			&boxsize,
			&cellsize,
			&(tot_x),&(tot_y)
			))  
    return NULL;

  //read from file(s) in a loop
  //***keep track of global particle index through pc and local file index through pc_new***
  pc = 0;
  pc_new = 0;
  countp = 0;

  for(i=0;i<nfiles;i++)
    {

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
	  //printf("allocate array N = %d\n",NumPart);
	  data = (float *) malloc(NumPart*3*sizeof(float));


	}
      
      for(k=0,ntot_withmasses=0;k<5;k++)
	if(header.mass[k] == 0)
	  ntot_withmasses += header.npart[k];

      fread(&dummy, sizeof(dummy), 1, fd);
      for(k=0;k<6;k++)
	{ 
	  //printf("Read particle %d : N = %d\n",k,header.npart[k]);
	  //printf("Read byte %d\n",fread(&(data[3*countp]),3*sizeof(float),header.npart[k],fd));
	  fread(&(data[3*countp]),3*sizeof(float),header.npart[k],fd);
	   
	  countp += header.npart[k];
	}
      fread(&dummy, sizeof(dummy), 1, fd);
      fclose(fd);
    }
  printf("Finish reading particles\n");
  radius = (xmax-xmin)/2.;
  //count = calloc(tot_x,sizeof(double *));
  countk = 0;
  for(i=0;i<tot_x;i++)
    {
      //count[i] = calloc(tot_y,sizeof(double));
    }
  dataout = malloc(2*(countk+1)*sizeof(float));
  for(i=0;i<NumPart;i++)
    {
      // Mpc => kpc
      // printf("%g\t%g\t%g\n",data[i*3],data[i*3+1],data[i*3+2]);
      rs[0] = data[i*3]*1000. - x0;
      rs[1] = data[i*3+1]*1000. - y0;
      rs[2] = data[i*3+2]*1000. -z0;
      for(j=0;j<3;j++)
	{
	  if(rs[j] > boxsize/2.)
	    rs[j] -= boxsize;
	  if(rs[j] < -1.*boxsize/2.)
	    rs[j] += boxsize;
	}
      length = rs[0]*z_axis[0]+rs[1]*z_axis[1]+rs[2]*z_axis[2];
      //printf("Length = %f\n",length);
      if(fabs(length) < radius)
	{
	  //printf("Length = %f\n",length);
	  for(j=0;j<3;j++)
	    {
	      rs[j] -= z_axis[j]*length;
	    }
	  x = rs[0]*x_axis[0]+rs[1]*x_axis[1]+rs[2]*x_axis[2];
	  y = rs[0]*y_axis[0]+rs[1]*y_axis[1]+rs[2]*y_axis[2];
	  if(x < xmax && x > xmin && y < ymax && y > ymin)
	    {
	      //printf("Length = %f\n",length);
	      //printf("%f %f\n",x,y);
	      
	      x_block = (int) ((x-xmin)/cellsize);
	      y_block = (int) ((y-ymin)/cellsize);
	      //printf("block %d %d\n",x_block, y_block);
	      //printf("size of float = %d\n",sizeof(float));
	      //printf("Going to allocate %ld elements = %ld bytes\n",2*(countk+1),2*(countk+1)*sizeof(float));
	      
	      dataout = realloc(dataout,2*(countk+1)*4);
	      
	      dataout[2*countk] = x;
	      dataout[2*countk+1] = y;
	      //printf("data = %g\t%g\n",dataout[2*(countk)],dataout[2*(countk)+1]);
	      //count[x_block][y_block]+= 1.;
	      countk++;
	      
	    }
	}

    }
  free(data);
  printf("finish transforming\n");
  //dim[0] = tot_x*tot_y;
  dim[0] = countk*2;
  rho = (float) NumPart/pow(boxsize,3)*(2.*radius*cellsize*cellsize);
  k=0;
  //dataout = calloc(tot_x*tot_y,sizeof(double));
  //smooth
  for(i=1;i<tot_x-1;i++)
    {
      for(j=1;j<tot_y-1;j++)
	{
	  //count[i][j] = count[i][j] + count[i+1][j] + count[i+1][j+1] + count[i+1][j-1] + count[i][j+1] + count[i][j-1] + count[i-1][j-1] +count[i-1][j] +count[i-1][j+1];
	  //count[i][j] /= 9.;
	}
    }
  k=0;
  for(i=0;i<tot_x;i++)
    {
      for(j=0;j<tot_y;j++)
	{
	  //dataout[k] =  (float) (count[i][j] / rho);
	  k++;
	}
    } 

  for(i=0;i<sizeof(dataout)/sizeof(dataout[0]); i+=2)
    {
      //printf("%f\t%f\n",dataout[i],dataout[i+1]);
    }
  printf("\n\n");
  printf("Allocate python object and copy data inz\n");
  //array = (PyArrayObject *) PyArray_FromDimsAndData(ndim, dim, NPY_DOUBLE,dataout);
  printf("total element %d %d\n",ndim,dim[0]/2);
  array = (PyArrayObject *)PyArray_SimpleNew(ndim, &dim,PyArray_FLOAT);
  //array = (PyArrayObject *) PyArray_Return(array);
  void *arr_data = PyArray_DATA((PyArrayObject*)array);
  memcpy(arr_data, dataout, PyArray_ITEMSIZE((PyArrayObject*) array) * dim[0]);
  return  (PyObject *)array;
}

PyMethodDef methods[] = {
  {"readone", gadget, METH_VARARGS, "Boyd"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC 
initgadgetPydensity()
{
  (void) Py_InitModule("gadgetPydensity", methods);   
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

